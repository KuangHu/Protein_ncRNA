#!/usr/bin/env python
"""Pfam domain annotation of the anchor proteins, per bag.

    pfam_annotate.py --domtbl <hmmscan domtblout> --out <dir> \\
                     [--nr100 bag_nr100.tsv] [--blocks blocks.tsv] \\
                     [--focus pilot_focus.txt]

This runs AFTER discovery and never feeds back into it. Discovery is anchored on
`cds_cluster_id` alone and knows nothing about what the protein is; this script
exists to say, in hindsight, what the anchors turned out to be. Nothing here is
in BLIND_CHAIN and nothing written here may be read by a discovery script -- if a
Pfam label ever selects which bags get carried forward, the blind is gone.

Domains are reported per *bag*, not per protein, because a bag is the unit the
pipeline carries: `n_reps_with_domain / n_reps` says how consistently the whole
bag is the thing the domain names, which is what distinguishes a bag that IS a
transposase from one that merely has a transposase representative.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

# Families worth calling out by name when they appear. These are read-outs, not
# filters: nothing is selected or dropped on the strength of this list, it only
# decides which column of the summary a hit is counted in. Keeping it explicit
# beats a substring search over descriptions, which would catch "transposase"
# inside unrelated prose.
GROUPS = {
    "IS110": ("Transposase_20", "Y1_Tnp", "DEDD_Tnp_IS110", "Transposase_IS110"),
    "group_II_intron": ("RVT_1", "RVT_2", "GIIM", "Intron_maturas2",
                        "LAGLIDADG_3", "Maturase"),
    "retron_RT": ("RVT_3", "RT_like"),
    "DDE_transposase": ("DDE_Tnp_1", "DDE_Tnp_1_2", "DDE_Tnp_1_3",
                        "DDE_Tnp_1_4", "DDE_Tnp_1_5", "DDE_Tnp_1_6",
                        "DDE_Tnp_1_7", "DDE_Tnp_IS1595", "DDE_Tnp_IS240",
                        "DDE_Tnp_ISL3", "DDE_3", "DDE_5"),
    "serine_recombinase": ("Resolvase", "Recombinase", "SpoIVCA"),
    "tyrosine_recombinase": ("Phage_integrase", "Phage_int_SAM_1", "Arm-DNA-bind"),
    "cas_nuclease": ("Cas9_REC", "HNH_Cas9", "RuvC_III", "Cas12", "TnpB_IS605"),
    "abortive_infection": ("Abi_2", "AbiEii", "AbiH", "Abi_HHH"),
}
OF = {d: g for g, ds in GROUPS.items() for d in ds}


def read_domtbl(p: Path):
    """Yield (query, target_name, target_acc, description, ievalue, score).

    hmmscan's domtblout is whitespace-delimited with the description running to
    end of line, so the split is capped at 22 fields rather than greedy.
    """
    with open(p) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split(None, 22)
            if len(f) < 23:
                continue
            yield f[3], f[0], f[1], f[22].strip(), float(f[12]), float(f[13])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domtbl", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--nr100", type=Path)
    ap.add_argument("--blocks", type=Path)
    ap.add_argument("--clade-stats", type=Path,
                    help="the step-4 clade_stats.tsv. Only bags listed there "
                         "were scored at all; without it, `n_real_blocks` "
                         "cannot distinguish 'no block' from 'never run'.")
    ap.add_argument("--focus", type=Path)
    ap.add_argument("--top", type=int, default=30)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    # Best hit per (protein, domain family); hmmscan reports one row per
    # envelope, and a tandem repeat of one domain is still one annotation.
    best: dict[tuple[str, str], tuple[float, float, str, str]] = {}
    prots_seen: set[str] = set()
    for q, name, acc, desc, iev, sc in read_domtbl(args.domtbl):
        prots_seen.add(q)
        k = (q, name)
        if k not in best or sc > best[k][1]:
            best[k] = (iev, sc, acc, desc)

    # The rep ids are "<bag>|<n>", so the bag is everything before the pipe.
    bag_of = lambda pid: pid.split("|")[0]
    by_bag_dom: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set))
    desc_of: dict[str, str] = {}
    acc_of: dict[str, str] = {}
    for (q, name), (iev, sc, acc, desc) in best.items():
        by_bag_dom[bag_of(q)][name].add(q)
        desc_of[name] = desc
        acc_of[name] = acc

    # Every rep in the FASTA, including the ones that matched nothing: a bag
    # whose reps are mostly un-annotated is a different finding from one whose
    # reps all hit the same domain, and only the denominator shows that.
    reps: dict[str, set[str]] = defaultdict(set)
    if args.nr100 and args.nr100.exists():
        pass
    for q in prots_seen:
        reps[bag_of(q)].add(q)
    # n_reps comes from the domtbl alone unless a rep list is supplied, so it is
    # "reps with at least one Pfam hit" -- stated plainly in the column name.

    nr100 = {}
    if args.nr100 and args.nr100.exists():
        with open(args.nr100) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                nr100[r["bag_id"]] = r

    # A bag that step 4 never saw must not be reported as a bag that produced no
    # block. The pilot scored 6 of 495 bags, so a bare `0` in that column would
    # read as a negative result for 489 bags that were never asked the question.
    # `scored` is the universe step 4 actually ran on; everything else is "".
    scored: set[str] = set()
    if args.clade_stats and args.clade_stats.exists():
        with open(args.clade_stats) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("decoy") == "False":
                    scored.add(r["clade"])
    has_block: dict[str, int] = defaultdict(int)
    if args.blocks and args.blocks.exists():
        with open(args.blocks) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("decoy") == "False":
                    has_block[r["clade"]] += 1
                    scored.add(r["clade"])

    focus = set(args.focus.read_text().split()) if args.focus else set()

    rows = []
    for bag in sorted(by_bag_dom):
        n = len(reps[bag])
        doms = by_bag_dom[bag]
        top = sorted(doms, key=lambda d: (-len(doms[d]), d))
        groups = sorted({OF[d] for d in doms if d in OF})
        rows.append({
            "bag": bag,
            "n_reps_with_a_pfam_hit": n,
            "n_distinct_domains": len(doms),
            "top_domain": top[0],
            "top_domain_acc": acc_of[top[0]],
            "top_domain_frac_of_reps": round(len(doms[top[0]]) / n, 4),
            "top_domain_desc": desc_of[top[0]],
            "all_domains": ";".join(
                f"{d}({len(doms[d])})" for d in top[:args.top]),
            "groups": ";".join(groups),
            "n_distinct_proteins": nr100.get(bag, {}).get("n_distinct", ""),
            "n_proteins": nr100.get(bag, {}).get("n_proteins", ""),
            "step4_scored": bag in scored,
            "n_real_blocks": has_block.get(bag, 0) if bag in scored else "",
            "in_pilot_focus": bag in focus,
        })

    cols = list(rows[0]) if rows else []
    out = args.out / "bag_pfam.tsv"
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # Catalogue-level tally: how much of the corpus each group accounts for.
    tally: dict[str, int] = defaultdict(int)
    for r in rows:
        for g in (r["groups"].split(";") if r["groups"] else ["(unclassified)"]):
            tally[g] += 1
    dom_tally: dict[str, int] = defaultdict(int)
    for bag, doms in by_bag_dom.items():
        for d in doms:
            dom_tally[d] += 1

    # What the ungrouped bags are, by their own top domain. Without this the
    # "(unclassified)" count is just a hole in the tally, and a reader cannot
    # tell a corpus of novel proteins from one whose families simply were not
    # on the GROUPS list.
    uncl: dict[str, int] = defaultdict(int)
    for r in rows:
        if not r["groups"]:
            uncl[r["top_domain"]] += 1

    (args.out / "bag_pfam_summary.json").write_text(json.dumps({
        "n_bags_with_a_hit": len(rows),
        "n_bags_scored_by_step4": sum(1 for r in rows if r["step4_scored"]),
        "n_proteins_with_a_hit": len(prots_seen),
        "top_domains_of_ungrouped_bags": dict(sorted(
            uncl.items(), key=lambda kv: -kv[1])[:40]),
        "bags_per_group": dict(sorted(tally.items(), key=lambda kv: -kv[1])),
        "bags_per_domain": dict(sorted(
            dom_tally.items(), key=lambda kv: -kv[1])[:60]),
    }, indent=2))

    print(f"{len(rows)} bags annotated, {len(prots_seen)} reps with a Pfam hit")
    print("\n== bags per group ==")
    for g, c in sorted(tally.items(), key=lambda kv: -kv[1]):
        print(f"  {c:5d}  {g}")
    print("\n== top domains by bag count ==")
    for d, c in sorted(dom_tally.items(), key=lambda kv: -kv[1])[:25]:
        print(f"  {c:5d}  {d:22s} {acc_of.get(d,''):12s} {desc_of.get(d,'')}")
    print("\n== top domains among ungrouped bags ==")
    for d, c in sorted(uncl.items(), key=lambda kv: -kv[1])[:20]:
        print(f"  {c:5d}  {d:22s} {desc_of.get(d,'')}")
    print(f"\n{sum(1 for r in rows if r['step4_scored'])} of {len(rows)} "
          f"annotated bags were scored by step 4; the rest have no "
          f"n_real_blocks value because they were never run.")
    if focus:
        print("\n== pilot focus bags ==")
        for r in rows:
            if r["in_pilot_focus"]:
                print(f"  {r['bag']}  blocks={r['n_real_blocks']}  "
                      f"{r['top_domain']} ({r['top_domain_frac_of_reps']}) "
                      f"[{r['groups'] or '-'}]  {r['top_domain_desc']}")
    print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
