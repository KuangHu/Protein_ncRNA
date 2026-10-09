#!/usr/bin/env python
"""Rank hard-null candidates on a continuous evidence scale.

    rank_candidates.py --candidates <discovery_candidates.tsv> \\
                       --null <null_comparison.tsv> --pool <pool_summary.tsv> \\
                       --floor 0.7077 --out <dir>

This changes no score. Every number it consumes was produced by the frozen
pipeline; this orders them and says how much each ordering rests on.

A pass/fail call on a null throws away most of what the run measured. "No decoy
produced a block" and "a decoy came within 0.06" are both nominally clean, and
they are not the same result; neither is "the null had six chances and found
nothing" the same as "no null could be built at all". So three quantities are
reported separately and never collapsed:

  signal      block quality, independent of any test. Recurrence, scaled by the
              number of DISTINCT anchor proteins it was measured over -- 0.99
              across 416 lineages and 0.99 across 34 are not equal evidence.
              Note that composite_score is very nearly recurrence alone for
              these bags: coord_tightness, strand_consistency and
              length_consistency all saturate at 1.0, so blending in
              top_member_fraction would double-count the same quantity.
  margin      how far the real block sits above the best challenge that was
              actually mounted -- the decoy maximum where decoys produced
              blocks, and otherwise the cross-arm floor. Which basis was used
              travels with the number.
  test        how much null pressure was really applied. A null that produced
              blocks is demonstrably capable of producing them; a null that
              produced none across six clades is weak evidence, not strong; no
              null at all is no evidence.

`evidence` is what the data demonstrates, `potential` is what the block would be
worth had it been fully tested and held up, and `gap` is the difference. A high
gap is not a bad candidate -- it is an untested one, and it is exactly where
further work buys the most. Ranking by evidence alone would bury the orphans,
which is the opposite of what the gap is for.

The three constants below are judgement calls about presentation, not pipeline
parameters, and are exposed as flags so a different weighting can be tried
without touching anything that scores.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def f(x, default=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def i(x, default=0):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return default


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True, type=Path)
    ap.add_argument("--null", required=True, type=Path)
    ap.add_argument("--pool", required=True, type=Path)
    ap.add_argument("--pfam", type=Path)
    ap.add_argument("--pfam-sensitive", type=Path,
                    help="pfam_sensitive.tsv; adds the sub-threshold status so "
                         "a true orphan is distinguishable from a bag that "
                         "merely has no domain above the curated bar")
    ap.add_argument("--targets", required=True, type=Path)
    ap.add_argument("--floor", type=float, required=True)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--null-mode", default="homology_pool_permutation")
    ap.add_argument("--breadth-ref", type=int, default=500,
                    help="distinct-anchor count treated as full lineage "
                         "breadth; below it the signal is discounted")
    ap.add_argument("--breadth-floor", type=float, default=0.5,
                    help="smallest breadth discount, so a narrow but real "
                         "block is reduced rather than erased")
    ap.add_argument("--margin-full", type=float, default=0.20,
                    help="separation above the best challenge that counts as a "
                         "full-credit margin")
    ap.add_argument("--untested-credit", type=float, default=0.25,
                    help="share of signal an untested block keeps; it is a "
                         "real recurrent block, so this is not zero")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    targets = args.targets.read_text().split()
    cand = {r["bag"]: r for r in csv.DictReader(
        open(args.candidates), delimiter="\t")}
    pool = {r["target"]: r for r in csv.DictReader(
        open(args.pool), delimiter="\t")}
    nullrow = {}
    for r in csv.DictReader(open(args.null), delimiter="\t"):
        if r["null_type"] == args.null_mode:
            nullrow[r["bag"]] = r
    doms = {}
    if args.pfam and args.pfam.exists():
        for r in csv.DictReader(open(args.pfam), delimiter="\t"):
            doms[r["bag"]] = r.get("all_domains", "")
    sens = {}
    if args.pfam_sensitive and args.pfam_sensitive.exists():
        sens = {r["bag"]: r for r in csv.DictReader(
            open(args.pfam_sensitive), delimiter="\t")}

    lo = math.log10(max(args.breadth_ref, 2))
    rows = []
    for t in targets:
        c, p = cand.get(t, {}), pool.get(t, {})
        n = nullrow.get(t, {})

        rec = f(c.get("max_composite_score"))
        n_anchors = i(c.get("n_distinct_anchors"))
        above = i(c.get("n_blocks_above_floor"))

        # Lineage breadth. log-scaled: the step from 30 to 300 distinct anchors
        # is worth far more than 300 to 3000.
        b = min(1.0, math.log10(max(n_anchors, 1)) / lo)
        breadth = args.breadth_floor + (1 - args.breadth_floor) * b
        signal = rec * breadth

        n_dc = i(n.get("n_decoy_clades"))
        n_db = i(n.get("n_decoy_blocks"))
        dmax = f(n.get("decoy_max_score"), float("nan"))

        if n_db > 0:
            challenge, basis = dmax, "decoy_max"
        else:
            challenge, basis = args.floor, "floor"
        margin = rec - challenge
        margin_factor = max(0.0, min(1.0, margin / args.margin_full))

        # A null that made blocks is proven able to make them. One that made
        # none may be a clean result or an inert null, and cannot tell us which.
        capability = 1.0 if n_db > 0 else (0.35 if n_dc > 0 else 0.0)
        test = capability * min(1.0, n_dc / 6.0)

        evidence = signal * margin_factor * (
            args.untested_credit + (1 - args.untested_credit) * test)
        potential = signal
        gap = potential - evidence

        if n_dc == 0:
            tier, tier_why = "D", "no null constructible from this corpus"
        elif n_db == 0:
            tier, tier_why = "C", "null built but produced no block of its own"
        elif margin >= 0.10:
            tier, tier_why = "A", "null demonstrably works and the block beat it"
        else:
            tier, tier_why = "B", "null demonstrably works and came close"

        rows.append({
            "rank": 0, "bag": t, "tier": tier,
            "evidence_score": round(evidence, 4),
            "potential_score": round(potential, 4),
            "information_gap": round(gap, 4),
            "recurrence": round(rec, 4),
            "n_distinct_anchors": n_anchors,
            "breadth_factor": round(breadth, 3),
            "signal": round(signal, 4),
            "n_blocks": i(c.get("n_blocks")),
            "n_blocks_above_floor": above,
            "n_decoy_clades": n_dc, "n_decoy_blocks": n_db,
            "decoy_max": "" if dmax != dmax else round(dmax, 4),
            "challenge_basis": basis,
            "challenge_level": round(challenge, 4),
            "margin": round(margin, 4),
            "margin_factor": round(margin_factor, 3),
            "test_strength": round(test, 3),
            "null_strength_source": p.get("null_strength_source", ""),
            "n_source_bags": i(p.get("n_source_bags")),
            "pool_distinct_proteins": i(p.get("pool_distinct_proteins")),
            "lane": c.get("lane", ""),
            "pfam_domains": doms.get(t, "") or "",
            "pfam_sensitive_status": sens.get(t, {}).get(
                "sensitive_status", ""),
            "pfam_sensitive_leads": sens.get(t, {}).get(
                "sensitive_leads", ""),
            "tier_reason": tier_why,
        })

    rows.sort(key=lambda r: (-r["evidence_score"], -r["signal"], r["bag"]))
    for k, r in enumerate(rows, 1):
        r["rank"] = k

    cols = list(rows[0]) if rows else []
    with open(args.out / "ranked_candidates.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    (args.out / "ranking_params.json").write_text(json.dumps({
        "floor": args.floor, "null_mode": args.null_mode,
        "breadth_ref": args.breadth_ref,
        "breadth_floor": args.breadth_floor,
        "margin_full": args.margin_full,
        "untested_credit": args.untested_credit,
        "note": "presentation weights; no pipeline score is affected",
    }, indent=2))

    print(f"{'#':>2s} {'bag':9s} {'T':1s} {'evid':>6s} {'pot':>6s} {'gap':>6s} "
          f"{'recur':>6s} {'anch':>5s} {'blk':>3s} {'dec':>3s} {'dblk':>4s} "
          f"{'margin':>7s} {'basis':10s} {'test':>5s}")
    for r in rows:
        print(f"{r['rank']:2d} {r['bag']:9s} {r['tier']:1s} "
              f"{r['evidence_score']:6.3f} {r['potential_score']:6.3f} "
              f"{r['information_gap']:6.3f} {r['recurrence']:6.4f} "
              f"{r['n_distinct_anchors']:5d} {r['n_blocks_above_floor']:3d} "
              f"{r['n_decoy_clades']:3d} {r['n_decoy_blocks']:4d} "
              f"{r['margin']:7.4f} {r['challenge_basis']:10s} "
              f"{r['test_strength']:5.2f}")
    for tier in "ABCD":
        n = sum(1 for r in rows if r["tier"] == tier)
        if n:
            why = next(r["tier_reason"] for r in rows if r["tier"] == tier)
            print(f"  tier {tier}: {n:2d}  {why}")
    top_gap = sorted(rows, key=lambda r: -r["information_gap"])[:4]
    print("highest information gap (most to gain from a better null): "
          + ", ".join(f"{r['bag']} {r['information_gap']:.3f}"
                      for r in top_gap))
    print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
