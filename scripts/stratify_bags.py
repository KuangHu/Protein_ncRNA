#!/usr/bin/env python
"""Stratify the Arm 2 bags by Pfam, then draw a balanced run list.

    stratify_bags.py --pfam <bag_pfam.tsv> --nr100 <bag_nr100.tsv> \\
                     --families <bag_families.tsv> --out <dir> \\
                     [--quota discovery=40,unknown=20,host_like=15]

Pfam is used here for STRATIFICATION, not prioritisation. The distinction is the
whole design: a bag is assigned to a lane before it is run, and the lane decides
which table its result is reported in, never whether it is allowed to score.

Four lanes:

  known_rna_mge   IS110, group II intron, IS200/IS605 (Y1), TnpB/Cas12f. These
                  are insertion systems that are already known to carry a
                  structured RNA. They are SPIKE-INS -- the measurement is
                  recovery rate, calibration for the blind caller. They are
                  excluded from novelty ranking by construction.
  discovery       transposases and recombinases with no known RNA partner.
                  This is the actual question.
  unknown         no Pfam hit at gathering thresholds, or nothing beyond a bare
                  HTH/DUF. Highest novelty, lowest interpretability.
  host_like       chemotaxis receptors, ushers, surface and metabolic proteins.
                  Not hard negatives -- they are the background estimate: if
                  these score as high as the transposases, the caller is broken.

Precedence is deliberate and conservative. A bag carrying both a known-RNA
marker and a generic transposase domain goes to `known_rna_mge`, so it cannot
leak into the discovery table and be reported as novel.

The domain lists below are a read-out vocabulary, not a model of the RNA. No
RNA family, length or structure appears anywhere here, and nothing in this file
is read by a discovery script.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

# Known RNA-carrying insertion systems, by the protein that marks them.
KNOWN_RNA_MGE = {
    "IS110": ("DEDD_Tnp_IS110", "Transposase_20"),
    "group_II_intron": ("GIIM", "RVT_1", "RVT_2", "RVT_N", "Intron_maturas2"),
    "IS200_IS605": ("Y1_Tnp",),
    "TnpB_Cas12": ("Cas12f1-like_TNB", "TnpB_IS605", "OrfB_IS605",
                   "HTH_OrfB_IS605", "Cas9_REC", "HNH_Cas9"),
}
# Mobile-element proteins with no known RNA partner. This is the lane the
# experiment is actually about.
DISCOVERY = {
    "DDE_transposase": (
        "DDE_Tnp_1", "DDE_Tnp_1_2", "DDE_Tnp_1_3", "DDE_Tnp_1_4", "DDE_Tnp_1_5",
        "DDE_Tnp_1_6", "DDE_Tnp_1_7", "DDE_Tnp_IS1595", "DDE_Tnp_IS240",
        "DDE_Tnp_ISL3", "DDE_Tnp_IS66", "DDE_Tnp_IS66_C", "DDE_Tnp_IS1",
        "DDE_Tnp_Tn3", "DDE_3", "DDE_5", "Zn_ribbon_ISL3", "zf-IS66",
        "LZ_Tnp_IS66", "Zn_ribbon_InsA", "DUF772"),
    "rve_integrase": ("rve", "rve_2", "rve_3", "Integrase_H2C2"),
    "Mu_transposase": ("Mu-transpos_C_2", "Mu-transpos_C", "Transpos_IS481"),
    "MULE_mutator": ("MULE", "Transposase_mut"),
    "tyrosine_recombinase": ("Phage_integrase", "Phage_int_SAM_1",
                             "Arm-DNA-bind"),
    "serine_recombinase": ("Resolvase", "Recombinase", "SpoIVCA"),
}
# Domains too generic to assign a lane on their own. A bag whose only hits are
# these is `unknown`, not a transposase: an HTH says the protein binds DNA, not
# what it is.
GENERIC = {"HTH_21", "HTH_23", "HTH_28", "HTH_38", "HTH_31", "HTH_17",
           "HTH_Tnp_1", "Helix-turn-helix", "zf-H2C2", "DUF4096"}

LANE_OF = {}
for grp, doms in KNOWN_RNA_MGE.items():
    for d in doms:
        LANE_OF[d] = ("known_rna_mge", grp)
for grp, doms in DISCOVERY.items():
    for d in doms:
        LANE_OF.setdefault(d, ("discovery", grp))


def assign(domains: list[str]) -> tuple[str, str]:
    """Lane and sub-group for one bag's domain list, known-RNA first."""
    hits = [LANE_OF[d] for d in domains if d in LANE_OF]
    known = [h for h in hits if h[0] == "known_rna_mge"]
    if known:
        # A bag can carry two markers (GIIM+RVT_1 is one system, but IS110 next
        # to an intron RT is two). Name them all rather than picking one.
        return "known_rna_mge", "+".join(sorted({g for _, g in known}))
    disc = [h for h in hits if h[0] == "discovery"]
    if disc:
        return "discovery", "+".join(sorted({g for _, g in disc}))
    if not domains or all(d in GENERIC for d in domains):
        return "unknown", "no_specific_domain"
    return "host_like", domains[0]


def spread(bags: list[dict], n: int) -> list[str]:
    """Take `n` bags by diversity, covering both family regimes.

    Sorting on `n_distinct` alone would fill every lane with singletons, since
    the most diverse Arm 2 bags are singletons. Half the budget is reserved for
    bags with >=4 siblings so the `same_family_permutation` null is available
    for part of each lane -- the pilot showed it is the only null that produces
    a near-miss, and a lane with none of it has no hard check.
    """
    key = lambda b: -int(b["n_distinct"] or 0)
    sib = sorted([b for b in bags if int(b["family_size"] or 1) >= 4], key=key)
    sol = sorted([b for b in bags if int(b["family_size"] or 1) < 4], key=key)
    out, i, j = [], 0, 0
    while len(out) < n and (i < len(sib) or j < len(sol)):
        if i < len(sib) and (len(out) % 2 == 0 or j >= len(sol)):
            out.append(sib[i]); i += 1
        elif j < len(sol):
            out.append(sol[j]); j += 1
    return [b["bag"] for b in out]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pfam", required=True, type=Path)
    ap.add_argument("--nr100", required=True, type=Path)
    ap.add_argument("--families", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--quota", default="discovery=40,unknown=20,host_like=15",
                    help="per-lane cap; known_rna_mge is always taken whole")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    quota = dict(kv.split("=") for kv in args.quota.split(","))
    quota = {k: int(v) for k, v in quota.items()}

    # The universe is the bags Arm 2 actually carries, not every bag with a
    # representative: bag_reps.faa holds 495, the Arm 2 corpus is the 399 in
    # bag_nr100.tsv. Using the larger set would schedule bags with no anchors.
    univ = {r["bag_id"]: r for r in csv.DictReader(
        open(args.nr100), delimiter="\t")}
    fam = {r["bag_id"]: r for r in csv.DictReader(
        open(args.families), delimiter="\t")}
    pfam = {r["bag"]: r for r in csv.DictReader(
        open(args.pfam), delimiter="\t")}

    rows = []
    for bag, nr in sorted(univ.items()):
        p = pfam.get(bag)
        doms = ([d.split("(")[0] for d in p["all_domains"].split(";")]
                if p and p["all_domains"] else [])
        lane, grp = assign(doms)
        f = fam.get(bag, {})
        rows.append({
            "bag": bag, "lane": lane, "lane_group": grp,
            "top_domain": p["top_domain"] if p else "",
            "top_domain_frac_of_reps": p["top_domain_frac_of_reps"] if p else "",
            "n_pfam_domains": len(doms),
            "all_domains": p["all_domains"] if p else "",
            "n_distinct": nr["n_distinct"], "n_proteins": nr["n_proteins"],
            "family": f.get("family", ""),
            "family_size": f.get("family_size", "1"),
        })

    by_lane: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_lane[r["lane"]].append(r)

    selected: dict[str, list[str]] = {}
    selected["known_rna_mge"] = sorted(
        r["bag"] for r in by_lane["known_rna_mge"])
    for lane in ("discovery", "unknown", "host_like"):
        selected[lane] = sorted(spread(by_lane[lane], quota.get(lane, 0)))

    chosen = {b: lane for lane, bs in selected.items() for b in bs}
    for r in rows:
        r["selected"] = r["bag"] in chosen

    cols = list(rows[0])
    with open(args.out / "bag_strata.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["lane"], r["bag"])):
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(args.out / "run_bags.txt", "w") as fh:
        for lane in ("known_rna_mge", "discovery", "unknown", "host_like"):
            for b in selected[lane]:
                fh.write(b + "\n")
    with open(args.out / "run_bag_lanes.tsv", "w") as fh:
        fh.write("bag\tlane\tlane_group\tn_distinct\tfamily_size\n")
        for r in sorted(rows, key=lambda r: (r["lane"], r["bag"])):
            if r["selected"]:
                fh.write(f"{r['bag']}\t{r['lane']}\t{r['lane_group']}\t"
                         f"{r['n_distinct']}\t{r['family_size']}\n")

    tot = sum(len(v) for v in selected.values())
    summary = {
        "n_bags_in_corpus": len(rows),
        "available_per_lane": {k: len(v) for k, v in sorted(by_lane.items())},
        "selected_per_lane": {k: len(v) for k, v in sorted(selected.items())},
        "n_selected": tot,
        "known_rna_mge_groups": dict(sorted(
            ((g, sum(1 for r in by_lane["known_rna_mge"]
                     if r["lane_group"] == g))
             for g in {r["lane_group"] for r in by_lane["known_rna_mge"]}),
            key=lambda kv: -kv[1])),
    }
    (args.out / "bag_strata_summary.json").write_text(json.dumps(summary, indent=2))

    print(f"{len(rows)} bags in corpus")
    for lane in ("known_rna_mge", "discovery", "unknown", "host_like"):
        print(f"  {lane:16s} available {len(by_lane[lane]):4d}   "
              f"selected {len(selected[lane]):4d}")
    print(f"  {'TOTAL':16s} {' ':14s}   selected {tot:4d}")
    print("\n== known_rna_mge spike-in lane ==")
    for r in sorted(by_lane["known_rna_mge"],
                    key=lambda r: (r["lane_group"], -int(r["n_distinct"]))):
        print(f"  {r['bag']}  {r['lane_group']:28s} "
              f"distinct={r['n_distinct']:>6s}  fam={r['family_size']:>3s}  "
              f"{r['top_domain']}")
    print(f"\n-> {args.out / 'run_bags.txt'}  ({tot} bags)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
