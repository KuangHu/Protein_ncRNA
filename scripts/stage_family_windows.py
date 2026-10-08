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
    ap.add_argument("--windows", required=True, type=Path)
    ap.add_argument("--families", required=True, type=Path)
    ap.add_argument("--subjects", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--min-size", type=int, default=4)
    args = ap.parse_args()
    (args.out / "proteins").mkdir(parents=True, exist_ok=True)

    fam = {r["bag_id"]: r["family"] for r in csv.DictReader(
        open(args.families), delimiter="\t")}
    size = Counter(fam.values())
    subjects = set(args.subjects.read_text().split())
    keep = {fam[b] for b in subjects
            if b in fam and size[fam[b]] >= args.min_size}
    print(f"{len(keep)} families with >= {args.min_size} members contain a "
          f"subject bag")

    n, skipped = 0, 0
    for d in sorted(args.windows.iterdir()):
        if not d.is_dir():
            continue
        bag = d.name
        f = fam.get(bag)
        if f not in keep:
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
