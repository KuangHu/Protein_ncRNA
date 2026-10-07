#!/usr/bin/env python
"""Draw a species-stratified genome sample out of genomes_db.

genomes_db holds 305k assemblies but only 9 species, and two of them are 79% of
the mirror. An unstratified sample is therefore almost entirely K. pneumoniae
and S. enterica. This takes an equal quota per species so that step 1-2 output
is at least comparable across hosts; it does NOT fix the underlying diversity
problem (see docs/PROJECT_STATUS.md).

It also skips **empty download stubs**, which are 84% of the mirror's view
entries (20-byte gzip files left by failed downloads). Sampling without this
filter wastes most of a job on genomes that make prodigal exit 18.
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

DB = Path("/global/scratch/users/kh36969/genomes_db")
VIEW = DB / "views/top10_IS110_total/by_species"
MIN_BYTES = 100_000  # a real bacterial assembly, gzipped, is well above this


def real_accessions(min_bytes: int) -> set[str]:
    """Accessions whose mirrored FASTA is a real assembly, per INDEX.tsv."""
    good = set()
    with open(DB / "INDEX.tsv") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[2].isdigit() and int(f[2]) >= min_bytes:
                good.add(f[0])
    return good


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-species", type=int, default=30)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min-bytes", type=int, default=MIN_BYTES)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    good = real_accessions(args.min_bytes)
    rng = random.Random(args.seed)
    rows, skipped = [], 0
    for sp_dir in sorted(VIEW.iterdir()):
        if not sp_dir.is_dir():
            continue
        files, stubs = [], 0
        for p in sorted(sp_dir.iterdir()):
            if not p.name.endswith(".fna.gz"):
                continue
            if p.name.replace(".fna.gz", "") in good:
                files.append(p)
            else:
                stubs += 1
        skipped += stubs
        pick = files if len(files) <= args.per_species else rng.sample(files, args.per_species)
        print(f"  {sp_dir.name:<30} {len(files):>6} real  {stubs:>6} stub  -> {len(pick)}")
        for p in sorted(pick):
            rows.append((p.name.replace(".fna.gz", ""), str(p.resolve()), sp_dir.name))

    with open(args.out, "w") as fh:
        fh.write("accession\tpath\tspecies\n")
        for acc, path, sp in rows:
            fh.write(f"{acc}\t{path}\t{sp}\n")
    print(f"\n{len(rows)} genomes from {len(set(r[2] for r in rows))} species "
          f"({skipped} stubs skipped) -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
