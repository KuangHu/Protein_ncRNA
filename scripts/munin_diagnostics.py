#!/usr/bin/env python
"""What is the RNAMunin score actually tracking?

A model that separates the RNA-carrying families only because their non-coding
regions are longer, or only because they are AT-rich, has not found an RNA. So
alongside the raw score this reports its relationship to length and GC, and
the AUC within length-matched strata.
"""
import argparse
import glob
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from separation_bench import POS, NEG, families, auc, perm_p, median, num


def read_fasta_lens(path):
    gc, ln, name = {}, {}, None
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                if name:
                    s = "".join(seq)
                    ln[name] = len(s)
                    gc[name] = (s.count("G") + s.count("C") +
                                s.count("g") + s.count("c")) / max(1, len(s))
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line)
    if name:
        s = "".join(seq)
        ln[name] = len(s)
        gc[name] = (s.count("G") + s.count("C") + s.count("g") +
                    s.count("c")) / max(1, len(s))
    return ln, gc


def spearman(x, y):
    def rank(v):
        o = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(o):
            j = i
            while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
                j += 1
            m = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[o[k]] = m
            i = j + 1
        return r
    rx, ry = rank(x), rank(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx) ** 0.5
    syy = sum((b - my) ** 2 for b in ry) ** 0.5
    return sxy / (sxx * syy) if sxx and syy else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root",
                    default="/global/scratch/users/kh36969/protein_ncrna/expand")
    a = ap.parse_args()

    lab = families(os.path.join(a.root, "pfam", "bag_pfam.tsv"))
    ln, gc = read_fasta_lens(os.path.join(a.root, "rnamunin",
                                          "nc_candidates.fasta"))

    rows = []
    with open(os.path.join(a.root, "rnamunin", "nc_munin_scores.tsv")) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            rows.append(dict(zip(head, line.rstrip("\n").split("\t"))))

    cand = {r["bag"]: r for r in rows if r["role"] == "candidate"}
    per = {b: max(num(r["mean_fwd"]), num(r["mean_rc"]))
           for b, r in cand.items()}
    lens = {b: float(r["len"]) for b, r in cand.items()}
    gcs = {b: gc.get(r["seq_id"], float("nan")) for b, r in cand.items()}

    have = [b for b in lab if b in per]
    print("labelled bags with a candidate score: %d" % len(have))
    xs = [per[b] for b in have]
    print("spearman(munin, length) = %.3f" % spearman(xs, [lens[b] for b in have]))
    print("spearman(munin, GC)     = %.3f" % spearman(xs, [gcs[b] for b in have]))

    pos = [per[b] for b in have if lab[b] in POS]
    neg = [per[b] for b in have if lab[b] in NEG]
    A, p = perm_p(pos, neg, 20000)
    print("munin candidate mean: AUC %.3f p=%.4f" % (A, p))
    print("GC alone:             AUC %.3f" %
          auc([gcs[b] for b in have if lab[b] in POS],
              [gcs[b] for b in have if lab[b] in NEG]))

    # length-matched strata: if the signal survives inside a length bin it is
    # not the length signal wearing a different hat
    order = sorted(have, key=lambda b: lens[b])
    q = len(order) // 4 or 1
    print("\nlength-stratified (quartiles of candidate length)")
    for i in range(4):
        grp = order[i * q:(i + 1) * q] if i < 3 else order[3 * q:]
        gp = [per[b] for b in grp if lab[b] in POS]
        gn = [per[b] for b in grp if lab[b] in NEG]
        if gp and gn:
            print("  len %4d-%-4d n+=%-3d n-=%-3d AUC %.3f"
                  % (lens[grp[0]], lens[grp[-1]], len(gp), len(gn),
                     auc(gp, gn)))
        else:
            print("  len %4d-%-4d n+=%-3d n-=%-3d --"
                  % (lens[grp[0]], lens[grp[-1]], len(gp), len(gn)))


if __name__ == "__main__":
    main()
