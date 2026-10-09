#!/usr/bin/env python
"""Summarise a sub-gathering-threshold Pfam scan as leads, not annotations.

    pfam_sensitive_report.py --domtbl <targets_sensitive.domtbl> \\
                             --strict <bag_pfam.tsv> --out <dir> \\
                             [--max-evalue 0.05] [--top 8]

The production scan uses --cut_ga, Pfam's curated per-family threshold, which is
the correct bar for saying a domain is present. This reads a scan run WITHOUT
that cutoff, where every hit is by construction below the curated bar.

The distinction it exists to preserve: for a bag with no curated domain, "no
domain above the gathering threshold" and "no detectable similarity to any Pfam
family" are entirely different statements, and only the second makes a protein
an orphan. Nothing here is an annotation, and the output column is named
`sensitive_leads` so it cannot be read as one.

Like pfam_annotate.py this is deliberately outside BLIND_CHAIN: it reads the
discovery out after the fact and must never select, rank or steer it.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domtbl", required=True, type=Path)
    ap.add_argument("--strict", required=True, type=Path)
    ap.add_argument("--targets", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--max-evalue", type=float, default=0.05,
                    help="i-Evalue above which a sub-threshold hit is treated "
                         "as noise and not listed as a lead")
    ap.add_argument("--top", type=int, default=8)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    strict = {r["bag"]: r.get("all_domains", "")
              for r in csv.DictReader(open(args.strict), delimiter="\t")}

    # domtblout is whitespace-separated with a free-text description last:
    # f[0] target, f[3] query, f[12] i-Evalue, f[13] domain score, f[22:] desc.
    best: dict[str, dict[str, tuple[float, float, str]]] = defaultdict(dict)
    n_reps: dict[str, set[str]] = defaultdict(set)
    for line in open(args.domtbl):
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 23:
            continue
        bag = f[3].split("|")[0]
        n_reps[bag].add(f[3])
        try:
            ie, sc = float(f[12]), float(f[13])
        except ValueError:
            continue
        cur = best[bag].get(f[0])
        if cur is None or ie < cur[0]:
            best[bag][f[0]] = (ie, sc, " ".join(f[22:]))

    bags = (args.targets.read_text().split() if args.targets
            else sorted(set(best) | set(strict)))
    rows = []
    for b in bags:
        hits = sorted(best.get(b, {}).items(), key=lambda kv: kv[1][0])
        keep = [h for h in hits if h[1][0] <= args.max_evalue][:args.top]
        has_strict = bool(strict.get(b, ""))
        # Orphan means nothing survives even without the curated bar. A bag with
        # curated domains is not an orphan whatever the sensitive scan shows.
        status = ("has_curated_domain" if has_strict else
                  "sub_threshold_lead" if keep else
                  "orphan_no_pfam_similarity")
        rows.append({
            "bag": b,
            "curated_domains": strict.get(b, ""),
            "sensitive_status": status,
            "n_sub_threshold_families": len(hits),
            "n_leads_at_evalue": len(keep),
            "best_lead": keep[0][0] if keep else "",
            "best_lead_ievalue": f"{keep[0][1][0]:.2g}" if keep else "",
            "best_lead_score": round(keep[0][1][1], 1) if keep else "",
            "best_lead_desc": keep[0][1][2] if keep else "",
            "sensitive_leads": ";".join(
                f"{n}({v[0]:.2g})" for n, v in keep),
        })

    cols = list(rows[0]) if rows else []
    with open(args.out / "pfam_sensitive.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    counts = defaultdict(int)
    for r in rows:
        counts[r["sensitive_status"]] += 1
    (args.out / "pfam_sensitive_summary.json").write_text(json.dumps({
        "max_evalue": args.max_evalue, "n_bags": len(rows),
        "status_counts": dict(counts),
        "orphans": [r["bag"] for r in rows
                    if r["sensitive_status"] == "orphan_no_pfam_similarity"],
    }, indent=2))

    print(f"{len(rows)} bags; " + ", ".join(
        f"{k}={v}" for k, v in sorted(counts.items())))
    for r in rows:
        if r["sensitive_status"] == "has_curated_domain":
            continue
        print(f"\n  {r['bag']}  {r['sensitive_status']}  "
              f"({r['n_sub_threshold_families']} families below the bar, "
              f"{r['n_leads_at_evalue']} at iE<={args.max_evalue})")
        for n, v in sorted(best.get(r["bag"], {}).items(),
                           key=lambda kv: kv[1][0])[:args.top]:
            if v[0] > args.max_evalue:
                continue
            print(f"      {n:22s} iE={v[0]:<9.2g} score={v[1]:6.1f}  "
                  f"{v[2][:58]}")
    print(f"\n-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
