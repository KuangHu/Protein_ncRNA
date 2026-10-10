#!/usr/bin/env python
"""Accept only blocks that lie inside the insertion, not in the host flank.

    insertion_filter.py --blocks <blocks.tsv> --block-fasta <blocks/> \\
                        --sites <bags_v4/insertions_sites.jsonl> \\
                        --out <dir> [--min-identity 80] [--min-cov 0.5]

The insertion is the mobile element and its architecture is `[non-coding][tnp]`
or `[tnp][non-coding]`. Only non-coding sequence INSIDE the element is a
candidate. The side is a property of each insertion, not a symmetry: where the
element is `[tnp][non-coding]`, everything upstream of the tnp is host genome
and must not be accepted at all.

The block caller does not know this. It carves a window +/-5 kb around the
anchor and searches both sides, so a conserved noncoding block it reports may
sit entirely in host flank -- recurrent, decoy-clean, above the floor, and still
not part of the mobile element. Score cannot separate the two; provenance can.

So each block is tested against the non-coding regions actually recorded inside
that bag's insertions (`noncoding_regions` in insertions_sites.jsonl). A block
that matches one is in the element. A block that matches none is host flank and
is rejected whatever it scored. This is an ACCEPTANCE criterion applied before
any ranking, not a feature added to one.

Both strands are searched on purpose: the `orient` field shows elements
inserting fwd and rc at the same site, so an element-internal region can appear
in either orientation relative to the window. That is orientation of the
element, not a free parameter.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path


def read_fasta(p: Path):
    name, buf, out = None, [], []
    for line in open(p):
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name:
        out.append((name, "".join(buf)))
    return out


def rep_seq(p: Path):
    fa = read_fasta(p)
    if not fa:
        return None
    fa.sort(key=lambda x: len(x[1]))
    return fa[len(fa) // 2][1]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--blocks", required=True, type=Path)
    ap.add_argument("--block-fasta", required=True, type=Path)
    ap.add_argument("--sites", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--pfam", type=Path)
    ap.add_argument("--min-identity", type=float, default=80.0)
    ap.add_argument("--min-cov", type=float, default=0.5,
                    help="aligned length as a fraction of the BLOCK, so a long "
                         "block is not accepted on a short shared motif")
    ap.add_argument("--max-nc-per-bag", type=int, default=3000)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    blocks = [r for r in csv.DictReader(open(args.blocks), delimiter="\t")
              if r["decoy"] == "False"]
    by_bag = defaultdict(list)
    for r in blocks:
        by_bag[r["clade"]].append(r)
    want = set(by_bag)
    print(f"{len(blocks)} real blocks across {len(want)} bags")

    # bag -> distinct in-insertion noncoding sequences
    nc = defaultdict(dict)
    nsites = defaultdict(int)
    with open(args.sites) as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            b = r.get("bag_id")
            if b not in want:
                continue
            nsites[b] += 1
            if len(nc[b]) >= args.max_nc_per_bag:
                continue
            for i, s in enumerate(r.get("noncoding_regions") or []):
                if len(s) >= 30:
                    nc[b].setdefault(s, f"{r.get('site_id','?')}|nc{i}")
    print(f"{sum(len(v) for v in nc.values())} distinct in-insertion noncoding "
          f"regions over {len(nc)} bags")
    missing = sorted(want - set(nc))
    if missing:
        print(f"WARNING {len(missing)} bags have no recorded insertion sites; "
              f"their blocks cannot be placed and are reported as "
              f"`unknown`, not rejected: {missing[:8]}")

    pf = {}
    if args.pfam and args.pfam.exists():
        pf = {r["bag"]: r.get("all_domains", "") for r in
              csv.DictReader(open(args.pfam), delimiter="\t")}

    tmp = Path(tempfile.mkdtemp())
    rows = []
    for bag in sorted(by_bag):
        regions = nc.get(bag, {})
        qs = {}
        for b in by_bag[bag]:
            s = rep_seq(args.block_fasta / f"{b['block_id']}.fasta")
            if s:
                qs[b["block_id"]] = s
        if not qs:
            continue
        best = {}
        if regions:
            q, sdb = tmp / "q.fa", tmp / "s.fa"
            q.write_text("".join(f">{k}\n{v}\n" for k, v in qs.items()))
            sdb.write_text("".join(
                f">nc{i}\n{s}\n" for i, s in enumerate(regions)))
            r = subprocess.run(
                ["blastn", "-task", "blastn", "-word_size", "11", "-dust", "no",
                 "-evalue", "1e-3", "-max_target_seqs", "5",
                 "-query", str(q), "-subject", str(sdb),
                 "-outfmt", "6 qseqid pident length bitscore"],
                capture_output=True, text=True)
            for line in r.stdout.splitlines():
                f = line.split("\t")
                if len(f) < 4:
                    continue
                qid, pid, alen, bits = f[0], float(f[1]), int(f[2]), float(f[3])
                if qid not in best or bits > best[qid][2]:
                    best[qid] = (pid, alen, bits)
        for b in by_bag[bag]:
            bid = b["block_id"]
            s = qs.get(bid)
            if s is None:
                continue
            pid, alen, _ = best.get(bid, (0.0, 0, 0.0))
            cov = alen / len(s) if s else 0.0
            if not regions:
                verdict = "unknown_no_sites"
            elif pid >= args.min_identity and cov >= args.min_cov:
                verdict = "in_insertion"
            else:
                verdict = "host_flank"
            rows.append({
                "block_id": bid, "bag": bag,
                "pfam_domains": pf.get(bag, ""),
                "side": b["side"], "dist_to_anchor": b["dist_to_anchor"],
                "block_len": len(s),
                "composite_score": b["composite_score"],
                "n_insertion_sites": nsites.get(bag, 0),
                "n_nc_regions": len(regions),
                "best_pident": round(pid, 1),
                "aligned_len": alen,
                "coverage_of_block": round(cov, 3),
                "verdict": verdict,
            })

    cols = list(rows[0]) if rows else []
    with open(args.out / "insertion_filter.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    tot = defaultdict(int)
    for r in rows:
        tot[r["verdict"]] += 1
    print("\nblock verdicts:")
    for k, v in sorted(tot.items(), key=lambda kv: -kv[1]):
        print(f"  {k:18s} {v:5d}  {100 * v / len(rows):5.1f}%")

    print("\nside of accepted vs rejected blocks:")
    for verdict in ("in_insertion", "host_flank"):
        c = defaultdict(int)
        for r in rows:
            if r["verdict"] == verdict:
                c[r["side"]] += 1
        print(f"  {verdict:14s} " + "  ".join(
            f"{k} {v}" for k, v in sorted(c.items())))

    bags_ok = sorted({r["bag"] for r in rows if r["verdict"] == "in_insertion"})
    (args.out / "bags_with_in_insertion_block.txt").write_text(
        "\n".join(bags_ok) + ("\n" if bags_ok else ""))
    print(f"\n{len(bags_ok)} of {len(by_bag)} bags keep at least one "
          f"in-insertion block")
    print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
