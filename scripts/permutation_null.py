#!/usr/bin/env python
"""Within-family permutation null: decoy clades that keep everything but the clade.

Blind — this never reads the benchmark branch.

    permutation_null.py --windows <windows_dir> --out <windows_dir>/decoys \\
                        --replicates 3 --seed 1

For each real selected clade of size N in family F, draw N windows at random
*without replacement from F's own pooled windows*, permuting the clade labels.
The decoy therefore matches the real clade in family, size, window length,
anchor geometry, ORF density, and genome/GC composition. The only thing taken
away is the one property step 4 is supposed to be detecting: that the members
are each other's relatives.

So any recurrent noncoding block step 4 reports in a decoy is an artifact — a
shared host gene, a mobile-element repeat, a composition bias, or a quirk of the
block caller. The decoy score distribution is the null that real blocks must
beat.

What this null does *not* control for is phylogenetic relatedness as such: decoy
members are unrelated, so they are a weak null in the sense that finding nothing
in them is unsurprising. It bounds the artifact rate, not the false-discovery
rate under realistic relatedness. A harder null (same clade, anchor coordinate
shuffled) is a separate control and is not this script.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_ncrna.seqio import read_fasta


def family_of(name: str) -> str:
    return name.split("__")[0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", required=True, type=Path,
                    help="directory of per-clade window FASTAs from extract_windows.py")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--replicates", type=int, default=3,
                    help="decoy clades drawn per real clade")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--min-size", type=int, default=8,
                    help="skip clades below the step-3 size floor")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    real = sorted(p for p in args.windows.glob("*.fasta"))
    if not real:
        raise SystemExit(f"no window FASTAs in {args.windows}")

    # Pool per family, keyed by anchor so the same window is never drawn twice
    # into one decoy even though it may sit in several real clades.
    pool: dict[str, dict[str, tuple[str, str]]] = defaultdict(dict)
    sizes: list[tuple[str, str, int]] = []
    for p in real:
        fam = family_of(p.stem)
        recs = [(h, s) for h, s in read_fasta(str(p))]
        for h, s in recs:
            pool[fam][h.split()[0]] = (h, s)
        if len(recs) >= args.min_size:
            sizes.append((fam, p.stem, len(recs)))

    for fam, members in sorted(pool.items()):
        print(f"{fam}: pool of {len(members)} distinct windows across "
              f"{sum(1 for f, _, _ in sizes if f == fam)} clades")

    manifest = []
    n_written = 0
    for fam, stem, n in sizes:
        keys = list(pool[fam])
        if len(keys) < n:
            print(f"  {stem}: pool {len(keys)} < clade size {n}, skipping",
                  file=sys.stderr)
            continue
        for rep in range(args.replicates):
            pick = rng.sample(keys, n)
            name = f"DECOY_{fam}__rep{rep}__{stem.split('__', 1)[1]}"[:180]
            with open(args.out / f"{name}.fasta", "w") as fh:
                for k in pick:
                    h, s = pool[fam][k]
                    fh.write(f">{h}\n{s}\n")
            manifest.append({"decoy": name, "family": fam, "mirrors": stem,
                             "replicate": rep, "n_windows": n})
            n_written += 1

    (args.out / "decoy_manifest.json").write_text(json.dumps(
        {"seed": args.seed, "replicates": args.replicates,
         "real_clades": len(sizes), "decoy_clades": n_written,
         "decoys": manifest}, indent=2))
    print(f"\n{n_written} decoy clades mirroring {len(sizes)} real clades -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
