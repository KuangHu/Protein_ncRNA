#!/usr/bin/env python
"""RNAMunin against every other IS family, not just the six-family negative.

The narrow benchmark (`separation_bench.py`) asked RNAMunin to rank IS110,
IS200/605 and group II above IS1, IS3 and IS4/5/10. It reached 0.656 and
IS110 alone came in at 0.660, which looked like a weak but real signal.

This repeats it against the negative `broad_is_bench.py` defines -- every bag
that looks like a transposase and carries no known non-coding RNA, twelve
families instead of three -- with the same aggregation sweep RNAMunin was
given before: per-record `mean` or `best50`, taken over a bag's records as
the mean, max or median, at minimum record lengths of 0, 100 and 200 bp.

Both strands are already in the score table; a record's value is the larger
of the two, as upstream scans both.

The per-family head-to-heads at the end are the point. An aggregate AUC over
a pooled negative can be carried by one easy family, and the question is
whether RNAMunin can tell an RNA-guided element from any particular kind of
RNA-free one.
"""
import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from separation_bench import auc, perm_p, median, num
from broad_is_bench import POS, label

RNA = ("IS110", "IS200/605", "GroupII")


def records(path):
    """bag -> [(length, mean, best50)], strand-max per record."""
    out = defaultdict(list)
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(head, line.rstrip("\n").split("\t")))
            m = max(num(f["mean_fwd"]), num(f["mean_rc"]))
            b = max(num(f["best50_fwd"]), num(f["best50_rc"]))
            out[f["bag"]].append((num(f["len"]) or 0, m, b))
    return out


AGG = {"mean": lambda v: sum(v) / len(v), "max": max, "med": median}


def aggregate(recs, metric, agg, minlen):
    i = 1 if metric == "mean" else 2
    out = {}
    for b, rs in recs.items():
        v = [r[i] for r in rs if r[0] >= minlen]
        if v:
            out[b] = AGG[agg](v)
    return out


def split(sc, lab, pos_families):
    pos = [v for b, v in sc.items()
           if lab.get(b) in pos_families]
    neg = [v for b, v in sc.items()
           if b in lab and lab[b] not in RNA]
    return pos, neg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root",
                    default="/global/scratch/users/kh36969/protein_ncrna/expand")
    ap.add_argument("--perm", type=int, default=20000)
    ap.add_argument("--cell",
                    help="pin the reported cell as metric,agg,minlen instead "
                         "of taking the sweep's best; the top two cells are "
                         "0.001 apart and the subsets they give differ")
    a = ap.parse_args()

    lab, skipped = label(os.path.join(a.root, "pfam", "bag_pfam.tsv"))
    recs = records(os.path.join(a.root, "rnamunin", "nc_munin_scores.tsv"))
    scored = {b for b in recs if b in lab}
    npos = sum(1 for b in scored if lab[b] in RNA)
    print("labelled %d bags; skipped %s" % (len(lab), dict(skipped)))
    print("RNAMunin records for %d labelled bags: %d RNA-guided, %d RNA-free"
          % (len(scored), npos, len(scored) - npos))
    print("records scored: %d" % sum(len(v) for b, v in recs.items()
                                     if b in lab))

    print("\n=== aggregation sweep, all three RNA-guided families ===")
    print("%-7s %-5s %-7s %5s %5s %7s" %
          ("metric", "agg", "minlen", "n+", "n-", "AUC"))
    best = None
    for metric in ("mean", "best50"):
        for agg in ("mean", "max", "med"):
            for minlen in (0, 100, 200):
                sc = aggregate(recs, metric, agg, minlen)
                pos, neg = split(sc, lab, RNA)
                v = auc(pos, neg)
                if v is None:
                    continue
                print("%-7s %-5s %-7d %5d %5d %7.3f"
                      % (metric, agg, minlen, len(pos), len(neg), v))
                if best is None or v > best[0]:
                    best = (v, metric, agg, minlen)

    v, metric, agg, minlen = best
    print("\nbest cell: %s / %s over a bag's records / minlen %d -> %.3f"
          % (metric, agg, minlen, v))
    if a.cell:
        metric, agg, minlen = a.cell.split(",")
        minlen = int(minlen)
        print("reporting the pinned cell instead: %s / %s / minlen %d"
              % (metric, agg, minlen))

    sc = aggregate(recs, metric, agg, minlen)
    print("\n=== at %s / %s / minlen %d ==="
          % (metric, agg, minlen))
    print("%-24s %4s %7s %9s" % ("RNA-guided set", "n", "AUC", "p"))
    sets = [("all three", RNA), ("without group II", ("IS110", "IS200/605"))]
    sets += [(f, (f,)) for f in RNA]
    for name, fams in sets:
        pos, neg = split(sc, lab, fams)
        o, p = perm_p(pos, neg, a.perm)
        print("%-24s %4d %7.3f %9.4f" % (name, len(pos), o, p))

    print("\n=== head to head, each RNA-guided family vs each RNA-free one ===")
    byfam = defaultdict(list)
    for b, s in sc.items():
        if b in lab:
            byfam[lab[b]].append(s)
    negf = sorted(f for f in byfam if f not in RNA)
    print("%-11s %s" % ("", " ".join("%9s" % f for f in negf)))
    for p in RNA:
        row = ["%9.3f" % auc(byfam[p], byfam[n]) if byfam[n] and byfam[p]
               else "%9s" % "-" for n in negf]
        print("%-11s %s" % (p, " ".join(row)))

    print("\n=== median score by family ===")
    for f, v in sorted(byfam.items(), key=lambda kv: -median(kv[1])):
        print("%-11s %6.3f  n=%d%s"
              % (f, median(v), len(v), "   <- RNA-guided" if f in RNA else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
