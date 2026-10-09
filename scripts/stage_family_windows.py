#!/usr/bin/env python
"""Restage per-bag windows under `<family>__<bag>` for permutation_null.py.

    stage_family_windows.py --windows <dir of per-bag dirs> \\
                            --families <bag_families.tsv> \\
                            --subjects <bag list> --out <dir> [--min-size 4]

permutation_null derives a clade's family from everything before the first `__`
in its filename, and draws each decoy from the *other* clades of that family.
Arm 2 windows are named `<bag>`, which would make every bag its own family with
a pool of one -- the degenerate null this staging exists to avoid. Copying them
as `<family>__<bag>` is what gives the sibling bags a pool to draw from.

Only families with at least `--min-size` members are staged, and only families
that contain a subject bag: staging a singleton reproduces the empty pool, and
staging a family no subject belongs to produces decoys nothing will use.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", required=True, type=Path, nargs="+",
                    help="one or more dirs of per-bag window dirs; searched in "
                         "order, first match wins")
    ap.add_argument("--families", type=Path,
                    help="bag_id -> family map (one family per bag)")
    ap.add_argument("--pairs", type=Path,
                    help="explicit family/bag_id pairs, for pools where a bag "
                         "legitimately belongs to several families")
    ap.add_argument("--subjects", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--min-size", type=int, default=4)
    args = ap.parse_args()
    if not (args.families or args.pairs):
        ap.error("one of --families or --pairs is required")
    (args.out / "proteins").mkdir(parents=True, exist_ok=True)

    # (family, bag) pairs are the general form; a flat map is the special case
    # where every bag has exactly one family.
    if args.pairs:
        pairs = [(r["family"], r["bag_id"]) for r in csv.DictReader(
            open(args.pairs), delimiter="\t")]
    else:
        pairs = [(r["family"], r["bag_id"]) for r in csv.DictReader(
            open(args.families), delimiter="\t")]
    size = Counter(f for f, _ in pairs)
    subjects = set(args.subjects.read_text().split())
    keep = {f for f, b in pairs
            if b in subjects and size[f] >= args.min_size}
    print(f"{len(keep)} families with >= {args.min_size} members contain a "
          f"subject bag")

    # One source dir per bag; later dirs are fallbacks for bags extracted in a
    # separate pass.
    src_of: dict[str, Path] = {}
    for root in args.windows:
        if not root.is_dir():
            continue
        for d in sorted(root.iterdir()):
            if d.is_dir() and (d / f"{d.name}.fasta").exists():
                src_of.setdefault(d.name, d)

    n, skipped = 0, 0
    for f, bag in sorted(set(pairs)):
        if f not in keep:
            continue
        d = src_of.get(bag)
        if d is None:
            print(f"  {bag}: no windows found, skipped", file=sys.stderr)
            skipped += 1
            continue
        jl, fa = d / f"{bag}.jsonl", d / f"{bag}.fasta"
        if not jl.exists() or not fa.exists():
            print(f"  {bag}: missing windows, skipped", file=sys.stderr)
            skipped += 1
            continue
        clade = f"{f}__{bag}"
        with open(jl) as src, open(args.out / f"{clade}.jsonl", "w") as dst:
            for line in src:
                r = json.loads(line)
                r["family"] = f
                dst.write(json.dumps(r) + "\n")
        with open(fa) as src, open(args.out / f"{clade}.fasta", "w") as dst:
            for line in src:
                if line.startswith(">"):
                    line = line.replace(f"family={bag}", f"family={f}", 1)
                dst.write(line)
        p = d / "proteins" / f"{bag}.faa"
        if p.exists():
            (args.out / "proteins" / f"{clade}.faa").write_bytes(p.read_bytes())
        n += 1
    staged = Counter(c.stem.split("__")[0] for c in args.out.glob("*.fasta"))
    thin = {f: c for f, c in staged.items() if c < args.min_size}
    print(f"{n} clades staged, {skipped} skipped -> {args.out}")
    if thin:
        # A family that lost members to extraction failures can fall below the
        # pool size permutation_null needs, and would then be dropped there
        # rather than here. Saying so now is the difference between a missing
        # null and an unexplained one.
        print(f"WARNING {len(thin)} staged families now hold fewer than "
              f"{args.min_size} clades and will yield no decoys: {thin}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
