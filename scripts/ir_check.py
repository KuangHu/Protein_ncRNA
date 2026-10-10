#!/usr/bin/env python
"""Is an anchor-tethered block a terminal inverted repeat rather than an RNA?

    ir_check.py --blocks <blocks.tsv> --block-fasta <blocks/> --pfam <tsv> \\
                --out <dir> [--min-identity 70] [--min-frac 0.5]

A DDE transposase is an insertion sequence, and an insertion sequence ends in
terminal inverted repeats: short, highly conserved, noncoding, at a fixed small
offset from the transposase ORF, present on BOTH sides. That is the same
description as a recurrent anchor-tethered noncoding block, and the block caller
has no way to tell them apart -- it was built to find conserved noncoding
segments tethered to a protein, and a TIR is exactly that. So "DDE transposase
bags carry blocks" is not by itself evidence of an RNA, and before any of them
is reported as one the alternative has to be tested rather than assumed.

Two tests, both on the called sequence itself:

  cross-strand   the best upstream block against the REVERSE COMPLEMENT of the
                 best downstream block of the same bag. A pair of terminal
                 inverted repeats matches here; two unrelated RNAs do not.
  self-palindrome   a block against its own reverse complement. A hairpin or a
                 single IR scores; a non-palindromic RNA does not.

A hit on either is not proof of a TIR, and a miss is not proof of an RNA -- a
structured RNA can be palindromic too, which is the point of a hairpin. What the
two separate the bags into is "has an alternative explanation that fits" and
"does not", and only the second group is worth carrying forward as a candidate.

Post-hoc: reads results, changes nothing, feeds nothing back.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

DDE_KEYS = ("DDE_Tnp", "DDE_3", "rve", "DUF772", "Tnp_ISL3", "DDE_Tnp_IS240",
            "DUF4372", "Mu-transpos")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


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


def best_hit(a: str, b: str, tmp: Path):
    """blastn a vs b on the PLUS strand only; best HSP (pident, alen, bits).

    -strand plus is not optional here. blastn searches both strands by default,
    so a query against its own reverse complement always returns a perfect
    full-length minus-strand hit -- it is the sequence aligned to itself. With
    both strands enabled this test called 96/96 DDE bags and 119/119 others
    palindromic, which is not a result, it is the default setting.
    """
    qa, qb = tmp / "a.fa", tmp / "b.fa"
    qa.write_text(f">a\n{a}\n")
    qb.write_text(f">b\n{b}\n")
    r = subprocess.run(
        ["blastn", "-task", "blastn", "-word_size", "7", "-dust", "no",
         "-strand", "plus", "-evalue", "1", "-query", str(qa),
         "-subject", str(qb), "-outfmt", "6 pident length bitscore"],
        capture_output=True, text=True)
    best = (0.0, 0, 0.0)
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) >= 3 and float(f[2]) > best[2]:
            best = (float(f[0]), int(f[1]), float(f[2]))
    return best


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--blocks", required=True, type=Path)
    ap.add_argument("--block-fasta", required=True, type=Path)
    ap.add_argument("--pfam", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--min-identity", type=float, default=70.0)
    ap.add_argument("--min-frac", type=float, default=0.5,
                    help="aligned length as a fraction of the shorter block")
    ap.add_argument("--only-in-insertion", type=Path,
                    help="insertion_filter.tsv; restrict to blocks the filter "
                         "placed inside the insertion. A host-flank block is "
                         "not a candidate at all, so asking whether it is an "
                         "inverted repeat answers a question nobody asked.")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    pf = {r["bag"]: r.get("all_domains", "") for r in
          csv.DictReader(open(args.pfam), delimiter="\t")}
    rows = [r for r in csv.DictReader(open(args.blocks), delimiter="\t")
            if r["decoy"] == "False"]
    if args.only_in_insertion and args.only_in_insertion.exists():
        keep = {r["block_id"] for r in csv.DictReader(
            open(args.only_in_insertion), delimiter="\t")
            if r["verdict"] == "in_insertion"}
        rows = [r for r in rows if r["block_id"] in keep]
        print(f"{len(rows)} blocks inside an insertion")
    by_bag = defaultdict(list)
    for r in rows:
        by_bag[r["clade"]].append(r)

    def seq_of(bid):
        p = args.block_fasta / f"{bid}.fasta"
        if not p.exists():
            return None
        fa = read_fasta(p)
        if not fa:
            return None
        fa.sort(key=lambda x: len(x[1]))
        return fa[len(fa) // 2][1]        # median-length member

    out = []
    tmpd = tempfile.mkdtemp()
    tmp = Path(tmpd)
    for bag, bs in sorted(by_bag.items()):
        dom = pf.get(bag, "")
        fam = "DDE/rve" if any(k in dom for k in DDE_KEYS) else "other"
        up = [b for b in bs if b["side"] == "upstream"]
        dn = [b for b in bs if b["side"] == "downstream"]
        up.sort(key=lambda b: -float(b["composite_score"]))
        dn.sort(key=lambda b: -float(b["composite_score"]))
        best = max(bs, key=lambda b: float(b["composite_score"]))
        sb = seq_of(best["block_id"])
        if sb is None:
            continue

        # self-palindrome: the block against its own reverse complement
        pid_s, alen_s, _ = best_hit(sb, rc(sb), tmp)
        pal = (pid_s >= args.min_identity
               and alen_s >= args.min_frac * len(sb))

        pid_x = alen_x = 0.0
        cross = False
        pair = ""
        if up and dn:
            su_, sd = seq_of(up[0]["block_id"]), seq_of(dn[0]["block_id"])
            if su_ and sd:
                pid_x, alen_x, _ = best_hit(su_, rc(sd), tmp)
                cross = (pid_x >= args.min_identity
                         and alen_x >= args.min_frac * min(len(su_), len(sd)))
                pair = f"{up[0]['block_id']}|{dn[0]['block_id']}"

        out.append({
            "bag": bag, "domain_class": fam,
            "pfam_domains": dom,
            "n_blocks": len(bs),
            "has_both_sides": bool(up and dn),
            "best_block": best["block_id"],
            "best_score": best["composite_score"],
            "best_len": len(sb),
            "dist_to_anchor": best["dist_to_anchor"],
            "self_palindrome": pal,
            "self_pident": round(pid_s, 1), "self_alen": alen_s,
            "cross_strand_pair": pair,
            "cross_strand_ir": cross,
            "cross_pident": round(pid_x, 1), "cross_alen": alen_x,
            "ir_explained": pal or cross,
        })

    cols = list(out[0]) if out else []
    with open(args.out / "ir_check.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"{'class':10s} {'bags':>5s} {'selfPal':>8s} {'crossIR':>8s} "
          f"{'either':>7s} {'either%':>8s}")
    for fam in ("DDE/rve", "other"):
        rs = [r for r in out if r["domain_class"] == fam]
        if not rs:
            continue
        sp = sum(1 for r in rs if r["self_palindrome"])
        cx = sum(1 for r in rs if r["cross_strand_ir"])
        ei = sum(1 for r in rs if r["ir_explained"])
        print(f"{fam:10s} {len(rs):5d} {sp:8d} {cx:8d} {ei:7d} "
              f"{100 * ei / len(rs):7.1f}%")
    print(f"\n-> {args.out}/ir_check.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
