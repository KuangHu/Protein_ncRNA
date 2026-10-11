#!/usr/bin/env python
"""Does a score put IS110, group II introns and IS200/605 above IS1/IS3/IS4-5-10?

One number per bag, one label per bag, AUC with a label-shuffle p-value, and
the same AUC again per RNA-carrying family against the pooled RNA-negative
bags, because a score that only sees introns is not a score for this problem.
"""
import argparse
import glob
import os
import random
import sys
from collections import defaultdict

POS = {"IS110": {"DEDD_Tnp_IS110", "Transposase_20"},
       "IS200/605": {"Y1_Tnp", "OrfB_IS605", "HTH_OrfB_IS605",
                     "Cas12f1-like_TNB"},
       "GroupII": {"RVT_1", "GIIM", "Intron_maturas2", "RVT_2"}}
NEG = {"IS1": {"DDE_Tnp_IS1", "HTH_Tnp_IS1"},
       "IS3": {"rve", "rve_2", "rve_3"},
       "IS4/5/10": {"DDE_Tnp_1", "DDE_Tnp_1_2", "DDE_Tnp_1_3", "DDE_Tnp_1_4",
                    "DDE_Tnp_1_5", "DDE_Tnp_1_6", "Nterm_IS4",
                    "DDE_Tnp_IS1595"}}


def families(pfam_tsv):
    """A bag is labelled only if its domains hit exactly one family; a bag
    that looks like two things at once is not evidence either way."""
    lab = {}
    with open(pfam_tsv) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ia, ib = head.index("all_domains"), head.index("bag")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            doms = {d.split("(")[0] for d in f[ia].split(";") if d}
            hit = [k for k, v in list(POS.items()) + list(NEG.items())
                   if doms & v]
            if len(hit) == 1:
                lab[f[ib]] = hit[0]
    return lab


def read_tsvs(pattern):
    rows = []
    for p in sorted(glob.glob(pattern)):
        with open(p) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            if not head or head == [""]:
                continue
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) == len(head):
                    rows.append(dict(zip(head, f)))
    return rows


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def auc(pos, neg):
    """Mann-Whitney U / (n*m) by midranks, so 20k shuffles stay cheap."""
    n, m = len(pos), len(neg)
    if not n or not m:
        return None
    allv = sorted((v, i >= n) for i, v in enumerate(pos + neg))
    rsum, i = 0.0, 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
            j += 1
        mid = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            if not allv[k][1]:
                rsum += mid
        i = j + 1
    return (rsum - n * (n + 1) / 2.0) / (n * m)


def perm_p(pos, neg, n=20000, seed=0):
    obs = auc(pos, neg)
    if obs is None:
        return None, None
    allv = pos + neg
    k = len(pos)
    rnd = random.Random(seed)
    ge = 0
    for _ in range(n):
        rnd.shuffle(allv)
        if auc(allv[:k], allv[k:]) >= obs:
            ge += 1
    return obs, (ge + 1.0) / (n + 1.0)


def median(v):
    v = sorted(v)
    if not v:
        return float("nan")
    m = len(v) // 2
    return v[m] if len(v) % 2 else (v[m - 1] + v[m]) / 2.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root",
                    default="/global/scratch/users/kh36969/protein_ncrna/expand")
    ap.add_argument("--perm", type=int, default=20000)
    ap.add_argument("--minerva", default="bench",
                    help="subdirectory of expand/minerva to read")
    a = ap.parse_args()

    lab = families(os.path.join(a.root, "pfam", "bag_pfam.tsv"))

    score = defaultdict(dict)

    for r in read_tsvs(os.path.join(a.root, "nc_insertion",
                                    "nc_candidates_*.tsv")):
        b = r["bag"]
        score[b]["cand_len"] = num(r.get("candidate_len"))
        score[b]["recurrence"] = num(r.get("recurrence"))
        ns, nr = num(r.get("n_sites")), num(r.get("n_nc_regions"))
        if ns:
            score[b]["nc_per_site"] = nr / ns

    mun = os.path.join(a.root, "rnamunin", "nc_munin_scores.tsv")
    if os.path.exists(mun):
        byb = defaultdict(lambda: defaultdict(list))
        with open(mun) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                f = dict(zip(head, line.rstrip("\n").split("\t")))
                b, role = f["bag"], f["role"]
                best = max(num(f["mean_fwd"]), num(f["mean_rc"]))
                b50 = max(num(f["best50_fwd"]), num(f["best50_rc"]))
                byb[b][role].append((best, b50))
                byb[b]["_all"].append((best, b50))
        for b, d in byb.items():
            if d.get("candidate"):
                score[b]["munin_cand_mean"] = d["candidate"][0][0]
                score[b]["munin_cand_best50"] = d["candidate"][0][1]
            allv = d["_all"]
            score[b]["munin_bag_mean"] = sum(x[0] for x in allv) / len(allv)
            score[b]["munin_bag_max50"] = max(x[1] for x in allv)
            score[b]["munin_bag_med"] = median([x[0] for x in allv])

    pat = os.path.join(a.root, "minerva", a.minerva, "minerva_regions_*.tsv")
    mv = defaultdict(list)
    for r in read_tsvs(pat):
        v = num(r.get("minerva_paired_fraction"))
        if v is not None:
            mv[r["bag"]].append(v)
    for b, v in mv.items():
        score[b]["minerva_paired"] = max(v)

    keys = ["cand_len", "recurrence", "nc_per_site", "munin_cand_mean",
            "munin_cand_best50", "munin_bag_mean", "munin_bag_med",
            "munin_bag_max50", "minerva_paired"]

    print("bags labelled: %d" % len(lab))
    for k in keys:
        have = [b for b in lab if score.get(b, {}).get(k) is not None]
        pos = [score[b][k] for b in have if lab[b] in POS]
        neg = [score[b][k] for b in have if lab[b] in NEG]
        if not pos or not neg:
            print("%-20s no data" % k)
            continue
        A, p = perm_p(pos, neg, a.perm)
        print("%-20s n+=%-4d n-=%-4d med+ %-9.4g med- %-9.4g AUC %.3f  p=%.4f"
              % (k, len(pos), len(neg), median(pos), median(neg), A, p))
        line = []
        for fam in POS:
            fp = [score[b][k] for b in have if lab[b] == fam]
            if fp:
                line.append("%s %d:%.3f" % (fam, len(fp), auc(fp, neg)))
        for fam in NEG:
            fn = [score[b][k] for b in have if lab[b] == fam]
            if fn:
                line.append("%s %d:med %.4g" % (fam, len(fn), median(fn)))
        print("    " + "  ".join(line))

    # a score that catches a family the other misses is worth multiplying
    combos = [("cand_len", "munin_bag_med"), ("cand_len", "munin_cand_mean"),
              ("cand_len", "minerva_paired"),
              ("munin_bag_med", "minerva_paired"),
              ("cand_len", "munin_bag_med", "minerva_paired")]
    print("\n--- rank-product combinations ---")
    for cmb in combos:
        have = [b for b in lab
                if all(score.get(b, {}).get(k) is not None for k in cmb)]
        if len(have) < 20:
            print("%-45s only %d bags" % ("x".join(cmb), len(have)))
            continue
        rank = {}
        for k in cmb:
            order = sorted(have, key=lambda b: score[b][k])
            for i, b in enumerate(order):
                rank[b] = rank.get(b, 0.0) + (i + 1.0) / len(order)
        pos = [rank[b] for b in have if lab[b] in POS]
        neg = [rank[b] for b in have if lab[b] in NEG]
        if not pos or not neg:
            print("%-45s no contrast" % "x".join(cmb))
            continue
        A, p = perm_p(pos, neg, a.perm)
        extra = []
        for fam in POS:
            fp = [rank[b] for b in have if lab[b] == fam]
            if fp:
                extra.append("%s %d:%.3f" % (fam, len(fp), auc(fp, neg)))
        print("%-45s n=%-4d n+=%-3d AUC %.3f p=%.4f  %s"
              % ("x".join(cmb), len(have), len(pos), A, p, "  ".join(extra)))


if __name__ == "__main__":
    main()
