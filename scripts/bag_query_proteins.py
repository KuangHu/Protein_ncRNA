#!/usr/bin/env python
"""Pick the representative proteins that Arm 2 searches the corpus with.

    bag_query_proteins.py --windows <bag windows dir> --out bag_reps.faa [-k 3]

One protein per bag would miss the bag's own diversity: a 500-member cluster at
30% identity spans more than any single member represents, and a query drawn
from one end of it recruits only that end. `k` representatives per bag,
`nearest_representatives`, so they spread across the bag without being its
outliers -- the farthest-point seeds themselves would be fragments and
misannotations, which is the opposite of what a query should be.

Ids are `<bag_id>|<n>`, which is what expand_by_homology.py splits on to decide
which bag a hit belongs to.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_ncrna.diversity import nearest_representatives
from protein_ncrna.seqio import read_fasta


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", required=True, type=Path,
                    help="bags_to_windows.py output (reads proteins/*.faa)")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("-k", "--per-bag", type=int, default=3)
    ap.add_argument("--min-aa", type=int, default=60)
    args = ap.parse_args()

    n_bags = n_out = 0
    with open(args.out, "w") as fh:
        for faa in sorted((args.windows / "proteins").glob("*.faa")):
            bag = faa.stem
            seqs = {h.split()[0]: s for h, s in read_fasta(str(faa))
                    if len(s) >= args.min_aa}
            if not seqs:
                continue
            n_bags += 1
            reps = nearest_representatives(sorted(seqs), seqs, args.per_bag)
            for i, r in enumerate(reps):
                fh.write(f">{bag}|{i}\n{seqs[r]}\n")
                n_out += 1
            print(f"  {bag}: {len(seqs)} proteins -> {len(reps)} reps", flush=True)
    print(f"\n{n_out} query proteins from {n_bags} bags -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
