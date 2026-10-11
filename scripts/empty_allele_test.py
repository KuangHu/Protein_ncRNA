#!/usr/bin/env python
"""Is a called block element-borne, or host sequence the element landed next to?

    empty_allele_test.py --blocks <blocks.tsv> --block-fasta <blocks/> \\
                         --empty <empty_alleles/> --index <index.tsv> \\
                         --out <dir> [--min-identity 90] [--min-cov 0.5]

The empty allele is the same locus in a genome that never received the element.
So the occupied/empty pair answers the question the window geometry could only
guess at:

    absent from the empty allele   -> the block arrived with the element
    present in the empty allele    -> the block is host sequence, and the
                                      element merely inserted beside it

This is the acceptance criterion from
`only-noncoding-inside-insertion`, but made from an observation instead of a
coordinate. `insertion_filter.py` asks whether a block matches a non-coding
region somebody already annotated, so it cannot tell novel element non-coding
from host flank -- both simply fail to match. This can, and it does not depend
on the element-extent estimate either, which is a per-bag median and is the
known source of the `beyond_extent` calls.

A block with no widened empty allele is reported `untested`, never rejected: the
absence of a reconstructed empty site is a fact about the genome sample, not
about the block.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def log(*a):
    print(*a, flush=True)


def require_slurm():
    if os.environ.get("PN_ALLOW_LOGIN_NODE") == "1":
        return
    if not os.environ.get("SLURM_JOB_ID"):
        sys.exit("refusing to run on a login node; submit with sbatch")


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
    ap.add_argument("--empty", required=True, type=Path)
    ap.add_argument("--index", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--min-identity", type=float, default=90.0)
    ap.add_argument("--min-cov", type=float, default=0.5,
                    help="aligned length as a fraction of the block")
    ap.add_argument("--max-alleles-per-bag", type=int, default=400)
    args = ap.parse_args()
    require_slurm()
    args.out.mkdir(parents=True, exist_ok=True)

    # bag -> the widened empty alleles of that bag's sites
    bag_of = {}
    for r in csv.DictReader(open(args.index), delimiter="\t"):
        if int(r["n_empty_alleles"] or 0) > 0:
            bag_of[r["site_id"]] = r["bag"]
    empties = defaultdict(list)
    for sid, bag in bag_of.items():
        p = args.empty / f"{sid}.fasta"
        if p.exists():
            empties[bag].append(p)
    log(f"{len(bag_of)} widened sites over {len(empties)} bags")

    blocks = [r for r in csv.DictReader(open(args.blocks), delimiter="\t")
              if r["decoy"] == "False"]
    by_bag = defaultdict(list)
    for r in blocks:
        by_bag[r["clade"]].append(r)

    rows = []
    tmp = args.out / "_tmp"
    tmp.mkdir(exist_ok=True)
    for bag in sorted(by_bag):
        qs = {}
        for b in by_bag[bag]:
            s = rep_seq(args.block_fasta / f"{b['block_id']}.fasta")
            if s:
                qs[b["block_id"]] = s
        if not qs:
            continue
        paths = empties.get(bag, [])[:args.max_alleles_per_bag]
        if not paths:
            for b in by_bag[bag]:
                if b["block_id"] in qs:
                    rows.append({"block_id": b["block_id"], "bag": bag,
                                 "n_empty_alleles": 0, "best_pident": 0.0,
                                 "coverage_of_block": 0.0,
                                 "verdict": "untested_no_empty_allele"})
            continue
        sdb = tmp / f"{bag}.empty.fa"
        n_alleles = 0
        with open(sdb, "w") as fh:
            for p in paths:
                for name, seq in read_fasta(p):
                    fh.write(f">{name}\n{seq}\n")
                    n_alleles += 1
        q = tmp / f"{bag}.blocks.fa"
        q.write_text("".join(f">{k}\n{v}\n" for k, v in qs.items()))
        res = subprocess.run(
            ["blastn", "-task", "blastn", "-word_size", "11", "-dust", "no",
             "-evalue", "1e-5", "-max_target_seqs", "50",
             "-query", str(q), "-subject", str(sdb),
             "-outfmt", "6 qseqid pident length bitscore"],
            capture_output=True, text=True).stdout
        best = {}
        for line in res.splitlines():
            f = line.split("\t")
            if len(f) < 4:
                continue
            qid, pid, alen, bits = f[0], float(f[1]), int(f[2]), float(f[3])
            if qid not in best or bits > best[qid][2]:
                best[qid] = (pid, alen, bits)
        for bid, s in qs.items():
            pid, alen, _ = best.get(bid, (0.0, 0, 0.0))
            cov = alen / len(s)
            present = pid >= args.min_identity and cov >= args.min_cov
            rows.append({
                "block_id": bid, "bag": bag, "n_empty_alleles": n_alleles,
                "block_len": len(s), "best_pident": round(pid, 1),
                "aligned_len": alen, "coverage_of_block": round(cov, 3),
                "verdict": "host_present_in_empty" if present
                           else "element_borne_absent_from_empty"})
        sdb.unlink(missing_ok=True)
        q.unlink(missing_ok=True)

    cols = ["block_id", "bag", "n_empty_alleles", "block_len", "best_pident",
            "aligned_len", "coverage_of_block", "verdict"]
    with open(args.out / "empty_allele_test.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    c = defaultdict(int)
    for r in rows:
        c[r["verdict"]] += 1
    log(f"\n{len(rows)} blocks tested against widened empty alleles:")
    for k, v in sorted(c.items(), key=lambda kv: -kv[1]):
        log(f"  {k:34s} {v:5d}  {100 * v / len(rows):5.1f}%")
    log(f"-> {args.out}/empty_allele_test.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
