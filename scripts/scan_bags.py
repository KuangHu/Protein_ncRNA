#!/usr/bin/env python
"""Open scan: rank every one-or-two-protein element by its non-coding RNA score.

Not a benchmark. Every bag whose element encodes one or two proteins and that
has element-trimmed windows is scored, labelled or not, by the two axes the
benchmark settled on -- how much conserved non-coding the element carries and
how much of it folds -- combined as a rank product.

The known RNA families ride along as positive controls. They are expected
near the top; a low rank for them would mean the scan is insensitive, not
that they are absent. The labels never enter the score.

Each row carries the two confounds that have already caught this score out:
`ir_arm`, the longest inverted repeat, because a transposon end pairs without
being transcribed, and `orf_z`, the longest stop-free stretch against
composition-matched shuffles, because unannotated protein is not an RNA.
"""
import argparse
import glob
import json
import os
import random
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from broad_is_bench import POS, NEG, label as broad_label
from ir_flag import longest_ir
from separation_bench import median, num

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
STOP = {"TAA", "TAG", "TGA"}


def longest_orf(s):
    best = 0
    for q in (s, s.translate(COMP)[::-1]):
        for f in range(3):
            run = 0
            for i in range(f, len(q) - 2, 3):
                if q[i:i + 3] in STOP:
                    best = max(best, run)
                    run = 0
                else:
                    run += 3
            best = max(best, run)
    return best


def orf_z(s, n=60, seed=0):
    s = s.upper()
    obs = longest_orf(s)
    r = random.Random(seed)
    c = list(s)
    v = []
    for _ in range(n):
        r.shuffle(c)
        v.append(longest_orf("".join(c)))
    m = sum(v) / len(v)
    sd = (sum((x - m) ** 2 for x in v) / len(v)) ** 0.5
    return (obs - m) / sd if sd > 1e-9 else 0.0


def read_candidates(cdir):
    out = {}
    for p in sorted(glob.glob(os.path.join(cdir, "*.fasta"))):
        bag = os.path.basename(p)[:-6]
        keep, buf = False, []
        with open(p) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    keep = "|candidate|" in line
                elif keep:
                    buf.append(line)
        if buf:
            out[bag] = "".join(buf)
    return out


def element_records(path):
    """median non-coding count, inserted length and coding length per bag"""
    agg = defaultdict(lambda: ([], [], []))
    with open(path) as fh:
        for line in fh:
            d = json.loads(line)
            a = agg[d["bag_id"]]
            for k, lst in (("nc_region_count", a[0]), ("inserted_len", a[1]),
                           ("nc_total_len", a[2])):
                try:
                    lst.append(int(d.get(k) or 0))
                except (TypeError, ValueError):
                    lst.append(0)
    return {b: (median(n), median(i), median(i) - median(t), len(n))
            for b, (n, i, t) in agg.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root",
                    default="/global/scratch/users/kh36969/protein_ncrna/expand")
    ap.add_argument("--sites",
                    default="/global/scratch/users/kh36969/fna_ins_discovery/"
                            "bags_v4/insertions_sites.jsonl")
    ap.add_argument("--max-cds", type=int, default=2)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    lab, _ = broad_label(os.path.join(a.root, "pfam", "bag_pfam.tsv"))
    cand = read_candidates(os.path.join(a.root, "nc_insertion", "candidates"))
    elem = element_records(a.sites)

    geo = {}
    with open(os.path.join(a.root, "element", "element_geometry.tsv")) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(h, line.rstrip("\n").split("\t")))
            geo[f["bag"]] = f

    # bench, broad and scan were run with identical settings, so they pool
    mv = defaultdict(list)
    for d in ("bench", "broad", "scan"):
        for p in glob.glob(os.path.join(a.root, "minerva", d,
                                        "minerva_regions_*.tsv")):
            with open(p) as fh:
                head = fh.readline().rstrip("\n").split("\t")
                for line in fh:
                    f = dict(zip(head, line.rstrip("\n").split("\t")))
                    v = num(f.get("minerva_paired_fraction"))
                    if v is not None:
                        mv[f["bag"]].append(v)

    rows = {}
    for bag, seq in cand.items():
        if bag not in elem or bag not in mv or bag not in geo:
            continue
        nrc, ins, coding, nsites = elem[bag]
        anchor = num(geo[bag].get("median_anchor_len")) or 0
        if anchor <= 0:
            continue
        extra = coding - anchor
        ncds = 1 if extra < 300 else (2 if extra < 1800 else 3)
        if ncds > a.max_cds:
            continue
        v = mv[bag]
        rows[bag] = {
            "n_cds": ncds, "cand_len": len(seq), "minerva": sum(v) / len(v),
            "n_regions": len(v), "elem_len": ins, "coding_len": coding,
            "nc_region_count": nrc, "n_sites": nsites,
            "ir_arm": longest_ir(seq) if len(seq) <= 600 else "",
            "orf_z": round(orf_z(seq), 2),
            "family": lab.get(bag, ""),
            "control": ("positive" if lab.get(bag) in POS else
                        ("negative" if lab.get(bag) in NEG else "")),
        }

    keys = list(rows)
    for k in ("cand_len", "minerva"):
        order = sorted(keys, key=lambda b: rows[b][k])
        for i, b in enumerate(order):
            rows[b].setdefault("rank", 0.0)
            rows[b]["rank"] += (i + 1.0) / len(order)
    for b in keys:
        rows[b]["score"] = rows[b]["rank"] / 2.0

    cols = ["bag", "score", "family", "control", "n_cds", "cand_len",
            "minerva", "ir_arm", "orf_z", "n_regions", "n_sites",
            "elem_len", "coding_len", "nc_region_count"]
    with open(a.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for b in sorted(keys, key=lambda x: -rows[x]["score"]):
            r = rows[b]
            r["bag"] = b
            out.write("\t".join(
                ("%.4f" % r[c]) if isinstance(r[c], float) else str(r[c])
                for c in cols) + "\n")
    sys.stderr.write("%d bags scanned -> %s\n" % (len(keys), a.out))


if __name__ == "__main__":
    main()
