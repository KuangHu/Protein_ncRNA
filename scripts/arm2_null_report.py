#!/usr/bin/env python
"""Compare Arm 2's candidate nulls on the same step-4 output.

    arm2_null_report.py --blocks <pilot/blocks> --decoys <pilot/decoys> \\
                        --focus <pilot_focus.txt> --out <dir>

One table per (bag, null type): how many decoy clades were built, how many
produced a block at all, where their scores land, and how far the real blocks
sit above them. The question is not which null gives the best-looking answer --
it is which one produces a null distribution at all, and whether the three
disagree about the same bag.

`null_strength` here is about RESOLVING POWER, not about how the decoys were
sourced. Arm 1 labelled every family `strong` on source-clade diversity while
1,497 decoy clades produced zero blocks between them, so every FDR cell was
vacuous and nothing in the label said so. A null that yields no blocks cannot
resolve any rate, and is reported as `no_resolution` however well-sourced it is.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

MODES = ["genomic_background", "within_bag_shuffle", "same_family_permutation"]


def read_tsv(p: Path) -> list[dict]:
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def pct(vals: list[float], f: float) -> float:
    if not vals:
        return float("nan")
    s = sorted(vals)
    return s[min(len(s) - 1, int(f * len(s)))]


def mode_of(clade: str) -> str | None:
    for m in MODES:
        if clade.startswith(f"DECOY_{m}__"):
            return m
    return None


def bag_of(clade: str, bags: set[str]) -> str | None:
    for b in bags:
        if clade.endswith(f"__{b}") or clade == b:
            return b
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--blocks", required=True, type=Path)
    ap.add_argument("--decoys", required=True, type=Path)
    ap.add_argument("--focus", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--floor", type=float, default=0.80)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    bags = set(args.focus.read_text().split())

    blocks = read_tsv(args.blocks / "blocks.tsv")
    S = lambda b: float(b["composite_score"])

    # Every decoy clade that was BUILT, not only those that produced a block.
    # The difference between the two is the whole measurement: a null whose
    # clades all came back empty looks identical to a perfect separation if only
    # the blocks are counted.
    built: dict[tuple[str, str], set[str]] = defaultdict(set)
    for p in sorted(args.decoys.glob("DECOY_*.fasta")):
        m = mode_of(p.stem)
        b = bag_of(p.stem, bags)
        if m and b:
            built[(b, m)].add(p.stem)

    real_scores: dict[str, list[float]] = defaultdict(list)
    decoy_scores: dict[tuple[str, str], list[float]] = defaultdict(list)
    decoy_with_block: dict[tuple[str, str], set[str]] = defaultdict(set)
    for blk in blocks:
        cl = blk["clade"]
        if blk["decoy"] == "True":
            m, b = mode_of(cl), bag_of(cl, bags)
            if m and b:
                decoy_scores[(b, m)].append(S(blk))
                decoy_with_block[(b, m)].add(cl)
        else:
            b = bag_of(cl, bags)
            if b:
                real_scores[b].append(S(blk))

    cols = ["bag", "null_type", "n_decoy_clades", "decoy_clades_with_a_block",
            "n_decoy_blocks", "decoy_max_score", "decoy_p90_score",
            "decoy_blocks_ge_floor", "n_real_blocks", "real_max_score",
            "real_blocks_ge_floor", "separation_real_max_minus_decoy_max",
            "empirical_fdr_ge_floor", "null_strength"]
    rows = []
    for b in sorted(bags):
        rs = real_scores.get(b, [])
        for m in MODES:
            nb = len(built[(b, m)])
            ds = decoy_scores[(b, m)]
            wb = len(decoy_with_block[(b, m)])
            r_ge = sum(1 for s in rs if s >= args.floor)
            d_ge = sum(1 for s in ds if s >= args.floor)
            if nb == 0:
                strength = "not_available"
            elif wb == 0:
                strength = "no_resolution"
            elif wb < 0.1 * nb:
                strength = "sparse"
            else:
                strength = "resolving"
            rows.append({
                "bag": b, "null_type": m, "n_decoy_clades": nb,
                "decoy_clades_with_a_block": wb, "n_decoy_blocks": len(ds),
                "decoy_max_score": round(max(ds), 4) if ds else "",
                "decoy_p90_score": round(pct(ds, 0.9), 4) if ds else "",
                "decoy_blocks_ge_floor": d_ge,
                "n_real_blocks": len(rs),
                "real_max_score": round(max(rs), 4) if rs else "",
                "real_blocks_ge_floor": r_ge,
                "separation_real_max_minus_decoy_max":
                    round(max(rs) - max(ds), 4) if rs and ds else "",
                # Suppressed unless the null actually resolved something, for
                # the Cas3/Arm 1 reason: a ratio against a null that produced
                # nothing is a missing measurement, not a rate of zero.
                "empirical_fdr_ge_floor":
                    round(d_ge / len(built[(b, m)]) / max(r_ge, 1), 4)
                    if r_ge and strength in ("resolving", "sparse") else "",
                "null_strength": strength,
            })

    out = args.out / "null_comparison.tsv"
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    (args.out / "null_comparison.json").write_text(json.dumps(
        {"floor": args.floor, "modes": MODES, "rows": rows}, indent=2))

    w = ["bag", "null_type", "n_decoy_clades", "decoy_clades_with_a_block",
         "decoy_max_score", "real_max_score",
         "separation_real_max_minus_decoy_max", "null_strength"]
    print("  ".join(c for c in w))
    for r in rows:
        print("  ".join(str(r[c]).ljust(len(c)) for c in w))
    print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
