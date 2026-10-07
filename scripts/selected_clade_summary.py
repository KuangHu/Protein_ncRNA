#!/usr/bin/env python
"""Per-family rollup of the selected clades — the step-4 entry-point summary.

    selected_clade_summary.py --clades <clades_dir>

Reads `selected_clades.tsv` and `clade_summary.tsv`, writes
`selected_clade_summary.tsv` and prints it. Kept separate from
clade_decompose.py so the rollup can be regenerated in seconds after a rule
change, without repeating the ~80 minute pairwise-identity pass.
"""

from __future__ import annotations

import argparse
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path

COLS = ["family", "threshold", "selected_clades", "complete_windows", "nr100",
        "median_identity", "median_aa_len", "species", "flagged_large",
        "flagged_single_species"]


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        cols = next(fh).rstrip("\n").split("\t")
        return [dict(zip(cols, line.rstrip("\n").split("\t"))) for line in fh]


def _f(v, default=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clades", required=True, type=Path)
    args = ap.parse_args()

    sel = read_tsv(args.clades / "selected_clades.tsv")
    if not sel:
        print("no selected clades")
        return 1

    by_fam: dict[tuple, list[dict]] = defaultdict(list)
    for r in sel:
        by_fam[(r["family"], r["threshold"])].append(r)

    rows = []
    for (fam, thr), rs in sorted(by_fam.items()):
        idents = [_f(r.get("median_identity")) for r in rs]
        idents = [i for i in idents if i is not None]
        species = {s for r in rs for s in [r.get("n_species", "0")]}
        rows.append({
            "family": fam,
            "threshold": thr,
            "selected_clades": len(rs),
            # The real size of the step-4 input: how many full +/-5 kb windows
            # are available to align, not how many anchors exist.
            "complete_windows": sum(int(r.get("n_complete_windows", 0) or 0) for r in rs),
            "nr100": sum(int(r.get("n_nr100", 0) or 0) for r in rs),
            "median_identity": round(statistics.median(idents), 4) if idents else "",
            "median_aa_len": int(statistics.median(
                [int(r.get("median_aa_len", 0) or 0) for r in rs])),
            "species": max(int(r.get("n_species", 0) or 0) for r in rs),
            "flagged_large": sum("large_clade" in (r.get("flags") or "") for r in rs),
            "flagged_single_species": sum(
                "single_species" in (r.get("flags") or "") for r in rs),
        })

    out = args.clades / "selected_clade_summary.tsv"
    with open(out, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in COLS) + "\n")

    w = [max(len(c), max((len(str(r[c])) for r in rows), default=0)) for c in COLS]
    print("  ".join(c.ljust(x) for c, x in zip(COLS, w)))
    for r in sorted(rows, key=lambda r: -r["complete_windows"]):
        print("  ".join(str(r[c]).ljust(x) for c, x in zip(COLS, w)))

    tot = {
        "families": len({r["family"] for r in rows}),
        "selected_clades": sum(r["selected_clades"] for r in rows),
        "complete_windows": sum(r["complete_windows"] for r in rows),
        "nr100": sum(r["nr100"] for r in rows),
    }
    print(f"\ntotal: {tot['selected_clades']} clades across {tot['families']} families, "
          f"{tot['complete_windows']} complete windows, {tot['nr100']} nonredundant anchors")
    (args.clades / "selected_clade_summary.json").write_text(
        json.dumps({"per_family": rows, "total": tot}, indent=2))
    print(f"-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
