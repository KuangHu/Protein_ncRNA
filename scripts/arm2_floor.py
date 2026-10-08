#!/usr/bin/env python
"""Recalibrate the score floor from Arm 2's own decoys.

    arm2_floor.py --blocks <blocks.tsv> --out <dir> [--lanes run_bag_lanes.tsv]

The 0.80 floor was set on the genomic corpus, with a different cap, a different
clade size distribution and a different meaning for `recurrence` -- this run
caps at each bag's distinct-protein count, so a member fraction is a fraction of
distinct anchors, not of raw copies. Carrying 0.80 across would be assuming the
answer. The floor here is whatever the decoys say it is.

Two floors are reported and neither is a p-value:

  floor_no_decoy_support   the highest composite score any decoy reached, over
                           the strict nulls only (within_bag_shuffle and
                           same_family_permutation). A real block above it has
                           NO OBSERVED DECOY SUPPORT. That is an empirical upper
                           bound set by the number of decoys actually scored --
                           it is not a false discovery rate, and with a few
                           hundred decoy clades it cannot be read as one.
  floor_p99                the 99th percentile of the strict-null decoy scores.
                           Less brittle to one outlying decoy, but it admits
                           decoys by construction, so anything between p99 and
                           the max is explicitly contaminated.

genomic_background is excluded from both. The pilot showed it produces no blocks
at all (0 of 18 bag-cells), so including it would drag any percentile toward
zero and manufacture a floor that nothing has to clear.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path

STRICT = ("within_bag_shuffle", "same_family_permutation")
MODES = ("genomic_background",) + STRICT


def mode_of(clade: str) -> str | None:
    for m in MODES:
        if clade.startswith(f"DECOY_{m}__"):
            return m
    return None


def pct(vals: list[float], f: float) -> float:
    if not vals:
        return float("nan")
    s = sorted(vals)
    return s[min(len(s) - 1, int(f * len(s)))]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--blocks", required=True, type=Path)
    ap.add_argument("--lanes", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    lane_of = {}
    if args.lanes and args.lanes.exists():
        lane_of = {r["bag"]: r["lane"] for r in csv.DictReader(
            open(args.lanes), delimiter="\t")}

    by_mode: dict[str, list[float]] = defaultdict(list)
    clades: dict[str, set[str]] = defaultdict(set)
    real: list[float] = []
    for r in csv.DictReader(open(args.blocks), delimiter="\t"):
        s = float(r["composite_score"])
        if r["decoy"] == "True":
            m = mode_of(r["clade"])
            if m:
                by_mode[m].append(s)
                clades[m].add(r["clade"])
        else:
            real.append(s)

    strict = [s for m in STRICT for s in by_mode[m]]
    n_strict_clades = sum(len(clades[m]) for m in STRICT)
    floor_max = max(strict) if strict else None
    floor_p99 = pct(strict, 0.99) if strict else None

    out = {
        "n_real_blocks": len(real),
        "n_strict_decoy_blocks": len(strict),
        "n_strict_decoy_clades_with_a_block": n_strict_clades,
        "floor_no_decoy_support": round(floor_max, 4) if strict else None,
        "floor_p99": round(floor_p99, 4) if strict else None,
        "legacy_genomic_floor": 0.80,
        "per_mode": {m: {
            "n_blocks": len(by_mode[m]),
            "n_clades_with_a_block": len(clades[m]),
            "max": round(max(by_mode[m]), 4) if by_mode[m] else None,
            "p99": round(pct(by_mode[m], 0.99), 4) if by_mode[m] else None,
            "p90": round(pct(by_mode[m], 0.90), 4) if by_mode[m] else None,
            "median": round(statistics.median(by_mode[m]), 4)
            if by_mode[m] else None,
        } for m in MODES},
    }
    if floor_max is not None:
        out["n_real_blocks_above_floor_no_decoy_support"] = sum(
            1 for s in real if s > floor_max)
        out["n_real_blocks_above_floor_p99"] = sum(
            1 for s in real if s > floor_p99)
    # A floor resting on a handful of decoy blocks is a weak floor, and the
    # number backing it travels with it rather than being looked up later.
    out["floor_backed_by_n_decoy_blocks"] = len(strict)
    if not strict:
        out["warning"] = (
            "no strict-null decoy produced a block; no floor can be calibrated "
            "from this run and nothing may be reported as decoy-clean")

    (args.out / "score_floor.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    print(f"\n-> {args.out / 'score_floor.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
