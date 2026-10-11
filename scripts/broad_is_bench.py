#!/usr/bin/env python
"""IS110 + IS200/605 + group II against every other IS family in the corpus.

The six-family benchmark asked a narrow question: three RNA carriers against
IS1, IS3 and IS4/5/10. This widens the negative to every bag that looks like
a transposase and has no known non-coding RNA -- ISL3, IS256, IS66, Mu, Tn3,
IS630 and the rest -- which is the population a blind caller actually runs
against.

Widening a negative class is where a benchmark gets quietly corrupted, so
anything that might carry an RNA is excluded rather than counted as a
negative: a reverse transcriptase domain (retrons, DGRs, group II), a Cas
domain, TnsB (Tn7-like elements guided by a crRNA), HNH (the HEARO RNAs) and
phage integrases, which are not IS elements at all.
"""
import argparse
import glob
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from separation_bench import auc, perm_p, median, num, read_tsvs

POS = {"IS110": {"DEDD_Tnp_IS110", "Transposase_20"},
       "IS200/605": {"Y1_Tnp", "OrfB_IS605", "HTH_OrfB_IS605",
                     "Cas12f1-like_TNB"},
       "GroupII": {"RVT_1", "GIIM", "Intron_maturas2", "RVT_2"}}

NEG = {"IS1": {"DDE_Tnp_IS1", "HTH_Tnp_IS1", "Zn_ribbon_InsA"},
       "IS3": {"rve", "rve_2", "rve_3"},
       "IS4/5/10": {"DDE_Tnp_1", "DDE_Tnp_1_2", "DDE_Tnp_1_3", "DDE_Tnp_1_4",
                    "DDE_Tnp_1_5", "DDE_Tnp_1_6", "Nterm_IS4",
                    "DDE_Tnp_IS1595", "DDE_Tnp_4", "DDE_Tnp_IS240"},
       "ISL3": {"DDE_Tnp_ISL3", "Zn_ribbon_ISL3", "HTH_Tnp_ISL3"},
       "IS256": {"Transposase_mut", "MULE"},
       "IS66": {"DDE_Tnp_IS66", "DDE_Tnp_IS66_C", "zf-IS66", "LZ_Tnp_IS66"},
       "Mu": {"Mu-transpos_C_2"},
       "Tn3": {"DDE_Tnp_Tn3", "DUF4158"},
       "IS630": {"DDE_3"}}

# a domain here means the bag is not usable as an RNA-free negative
UNSAFE = {"RVT_1", "RVT_2", "RVT_3", "GIIM", "Intron_maturas2",
          "TnsB_C", "TnsD", "TniQ", "HNH", "Phage_integrase",
          "DEDD_Tnp_IS110", "Transposase_20", "Y1_Tnp", "OrfB_IS605",
          "HTH_OrfB_IS605"}


def label(pfam_tsv):
    """A bag gets a family only if exactly one family's domains are present.
    Negatives additionally have to be free of every UNSAFE domain."""
    lab, skipped = {}, defaultdict(int)
    with open(pfam_tsv) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ia, ib = head.index("all_domains"), head.index("bag")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            doms = {d.split("(")[0] for d in f[ia].split(";") if d}
            hit = [k for k, v in list(POS.items()) + list(NEG.items())
                   if doms & v]
            if len(hit) != 1:
                skipped["ambiguous_or_none" if hit else "no_family"] += 1
                continue
            fam = hit[0]
            if fam in NEG and (doms & UNSAFE):
                skipped["negative_with_rna_capable_domain"] += 1
                continue
            if fam in POS and any(d in doms for d in ("Cas9", "Cas12k")):
                pass
            lab[f[ib]] = fam
    return lab, skipped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root",
                    default="/global/scratch/users/kh36969/protein_ncrna/expand")
    ap.add_argument("--perm", type=int, default=20000)
    ap.add_argument("--ir-cut", type=int, default=12)
    a = ap.parse_args()

    lab, skipped = label(os.path.join(a.root, "pfam", "bag_pfam.tsv"))
    print("labelled %d bags; skipped %s" % (len(lab), dict(skipped)))

    sc = defaultdict(dict)
    for r in read_tsvs(os.path.join(a.root, "nc_insertion",
                                    "nc_candidates_*.tsv")):
        sc[r["bag"]]["len"] = num(r["candidate_len"])
    # bench covers the six original families, broad the rest; identical
    # settings, so they pool
    mv = defaultdict(list)
    pats = [os.path.join(a.root, "minerva", d, "minerva_regions_*.tsv")
            for d in ("bench", "broad")]
    for p in sorted(q for pat in pats for q in glob.glob(pat)):
        with open(p) as fh:
            h = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                f = dict(zip(h, line.rstrip("\n").split("\t")))
                v = num(f.get("minerva_paired_fraction"))
                if v is not None:
                    mv[f["bag"]].append(v)
    for b, v in mv.items():
        sc[b]["minerva"] = sum(v) / len(v)
    ir = {}
    irf = os.path.join(a.root, "nc_insertion", "ir_flag.tsv")
    if os.path.exists(irf):
        with open(irf) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 5 and f[2]:
                    ir[f[0]] = int(f[2])

    have = [b for b in lab if "len" in sc[b] and "minerva" in sc[b]]
    print("\nbags with both scores: %d" % len(have))
    print("%-11s %-5s %-9s %-9s" % ("family", "n", "med len", "med minerva"))
    for fam in list(POS) + list(NEG):
        bs = [b for b in have if lab[b] == fam]
        if bs:
            print("%-11s %-5d %-9.0f %-9.4f"
                  % (fam, len(bs), median([sc[b]["len"] for b in bs]),
                     median([sc[b]["minerva"] for b in bs])))

    def run(sub, title):
        rank = {}
        for k in ("len", "minerva"):
            o = sorted(sub, key=lambda b: sc[b][k])
            for i, b in enumerate(o):
                rank[b] = rank.get(b, 0.0) + (i + 1.0) / len(o)
        pos = [rank[b] for b in sub if lab[b] in POS]
        neg = [rank[b] for b in sub if lab[b] in NEG]
        A, p = perm_p(pos, neg, a.perm)
        print("\n%s: n=%d  n+=%d  n-=%d  AUC %.3f  p=%.4f"
              % (title, len(sub), len(pos), len(neg), A, p))
        for fam in POS:
            fp = [rank[b] for b in sub if lab[b] == fam]
            if fp:
                print("   %-11s n=%-3d AUC %.3f" % (fam, len(fp), auc(fp, neg)))
        for fam in NEG:
            fn = [rank[b] for b in sub if lab[b] == fam]
            if fn:
                print("   %-11s n=%-3d median rank pct %.2f  above-pos-p10 %d"
                      % (fam, len(fn), median(fn) / 2,
                         sum(1 for x in fn
                             if x >= sorted(pos)[int(len(pos) * 0.10)])))
        return rank

    run(have, "all bags")
    if ir:
        sub = [b for b in have if ir.get(b, 99) < a.ir_cut]
        run(sub, "inverted-repeat arm < %d" % a.ir_cut)


if __name__ == "__main__":
    main()
