#!/usr/bin/env python
"""The three axes together: conserved length, folding, and strand consistency.

`broad_is_bench.py` scores a bag on how much conserved non-coding its element
carries (`len`) and how much of it folds (`minerva`). This adds the third
axis, which only became computable once `orient_candidates.py` put the
candidates into the transposase's frame: how consistently the candidate sits
on one strand of the protein across the windows it was found in.

A transcript has a strand. A conserved stretch of DNA that is not transcribed
does not have to, and a transposon end -- which pairs with its partner -- has
a reason to appear on both.

Rank sum over the axes, as in `broad_is_bench.py`, so the three rows of the
table are comparable. Every combination is evaluated on the SAME bags, the
ones carrying all three axes, because a score that silently drops bags can
gain AUC from the subset rather than from the feature.
"""
import argparse
import glob
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from separation_bench import auc, perm_p, median, num
from broad_is_bench import POS, NEG, label


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root",
                    default="/global/scratch/users/kh36969/protein_ncrna/expand")
    ap.add_argument("--perm", type=int, default=20000)
    ap.add_argument("--ir-cut", type=int, default=12)
    ap.add_argument("--min-hits", type=int, default=3,
                    help="window placements a strand fraction needs to mean "
                         "anything")
    a = ap.parse_args()

    lab, skipped = label(os.path.join(a.root, "pfam", "bag_pfam.tsv"))
    sc = defaultdict(dict)

    for p in sorted(glob.glob(os.path.join(a.root, "nc_insertion",
                                           "nc_candidates_*.tsv"))):
        with open(p) as fh:
            h = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                f = dict(zip(h, line.rstrip("\n").split("\t")))
                sc[f["bag"]]["len"] = num(f["candidate_len"])

    mv = defaultdict(list)
    for d in ("bench", "broad"):
        for p in sorted(glob.glob(os.path.join(a.root, "minerva", d,
                                               "minerva_regions_*.tsv"))):
            with open(p) as fh:
                h = fh.readline().rstrip("\n").split("\t")
                for line in fh:
                    f = dict(zip(h, line.rstrip("\n").split("\t")))
                    v = num(f.get("minerva_paired_fraction"))
                    if v is not None:
                        mv[f["bag"]].append(v)
    for b, v in mv.items():
        sc[b]["minerva"] = sum(v) / len(v)

    # strand consistency: how far the candidate's placements are from an even
    # split between the two strands of the transposase
    with open(os.path.join(a.root, "nc_insertion", "oriented",
                           "orientation.tsv")) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(h, line.rstrip("\n").split("\t")))
            hits, fp = num(f.get("n_window_hits")), num(f.get("frac_plus"))
            if fp is not None and (hits or 0) >= a.min_hits:
                sc[f["bag"]]["strand"] = max(fp, 1.0 - fp)

    ir = {}
    with open(os.path.join(a.root, "nc_insertion", "ir_flag.tsv")) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5 and f[2]:
                ir[f[0]] = int(f[2])

    axes = ("len", "minerva", "strand")
    have = [b for b in lab if all(k in sc[b] for k in axes)]
    pos_b = [b for b in have if lab[b] in POS]
    neg_b = [b for b in have if lab[b] in NEG]
    print("labelled %d bags; skipped %s" % (len(lab), dict(skipped)))
    print("bags with all three axes: %d  (%d RNA-guided, %d RNA-free)"
          % (len(have), len(pos_b), len(neg_b)))

    print("\n=== each axis alone, on those bags ===")
    for k in axes:
        o, p = perm_p([sc[b][k] for b in pos_b], [sc[b][k] for b in neg_b],
                      a.perm)
        print("%-9s AUC %.3f  p=%.4f   median +%.3f  -%.3f"
              % (k, o, p, median([sc[b][k] for b in pos_b]),
                 median([sc[b][k] for b in neg_b])))

    def run(sub, keys, title):
        rank = defaultdict(float)
        for k in keys:
            o = sorted(sub, key=lambda b: sc[b][k])
            for i, b in enumerate(o):
                rank[b] += (i + 1.0) / len(o)
        pos = [rank[b] for b in sub if lab[b] in POS]
        neg = [rank[b] for b in sub if lab[b] in NEG]
        A, p = perm_p(pos, neg, a.perm)
        per = []
        for fam in POS:
            fp = [rank[b] for b in sub if lab[b] == fam]
            per.append("%s %.3f" % (fam, auc(fp, neg)) if fp
                       else "%s -" % fam)
        print("%-38s n=%-4d AUC %.3f  p=%.4f   %s"
              % (title, len(sub), A, p, "  ".join(per)))

    clean = [b for b in have if ir.get(b, 99) < a.ir_cut]
    print("\n=== combined, identical bag sets within each block ===")
    run(have, ("len", "minerva"), "len x minerva")
    run(have, axes, "len x minerva x strand")
    print()
    run(clean, ("len", "minerva"), "len x minerva, IR arm < %d" % a.ir_cut)
    run(clean, axes, "len x minerva x strand, IR arm < %d" % a.ir_cut)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
