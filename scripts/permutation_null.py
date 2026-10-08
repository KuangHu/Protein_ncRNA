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

Drawing uniformly from the family pool is not enough to break relatedness, and
the failure is quiet. The pool is dominated by its largest clades, so a decoy
can be mostly one real clade wearing a decoy name -- measured on the 43-clade
set, Cas3 decoys were 67-71% drawn from a single real clade, and the median
decoy over all families drew 33% of its members from one source. A family with
one selected clade degenerates completely: its decoy is a permutation of itself.

So each decoy is built with a per-source quota: members are drawn round-robin
from the family's *other* clades, at most ceil(N / n_sources) from any one of
them, and the clade being mirrored is excluded outright.

The bias this removes runs in the conservative direction, which is why it was
invisible. A decoy that retains relatedness behaves like a real clade, produces
more blocks, and *raises* the measured FDR. A contaminated null therefore
overstates the false-discovery rate rather than understating it -- so a family
reporting zero decoy blocks is safe, while a family reporting many may simply
have had a degenerate null. Each decoy records `n_source_clades` and
`max_source_fraction`, and a family whose decoys cannot meet the quota is marked
`null_strength: degenerate` so its FDR is not quoted as if it were measured.
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
    ap.add_argument("--min-source-clades", type=int, default=4,
                    help="a decoy should mix at least this many source clades; "
                         "below it the null is marked degenerate")
    ap.add_argument("--max-source-fraction", type=float, default=0.5,
                    help="a decoy drawing more than this fraction from one "
                         "source clade is marked degenerate")
    ap.add_argument("--allow-self-source", action="store_true",
                    help="let a decoy draw from the clade it mirrors; off by "
                         "default because that is the most direct leak")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    real = sorted(p for p in args.windows.glob("*.fasta"))
    if not real:
        raise SystemExit(f"no window FASTAs in {args.windows}")

    # Pool per family, keyed by anchor so the same window is never drawn twice
    # into one decoy even though it may sit in several real clades, and tracked
    # back to its source clade so relatedness can be budgeted rather than hoped
    # away.
    pool: dict[str, dict[str, tuple[str, str]]] = defaultdict(dict)
    by_source: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    sizes: list[tuple[str, str, int]] = []
    for p in real:
        fam = family_of(p.stem)
        recs = [(h, s) for h, s in read_fasta(str(p))]
        for h, s in recs:
            k = h.split()[0]
            pool[fam][k] = (h, s)
            by_source[fam][p.stem].append(k)
        if len(recs) >= args.min_size:
            sizes.append((fam, p.stem, len(recs)))

    for fam, members in sorted(pool.items()):
        print(f"{fam}: pool of {len(members)} distinct windows across "
              f"{sum(1 for f, _, _ in sizes if f == fam)} clades")

    def draw(fam: str, stem: str, n: int) -> tuple[list[str], dict]:
        """n windows from fam, spread across its other clades under a quota."""
        srcs = {s: list(v) for s, v in by_source[fam].items()
                if args.allow_self_source or s != stem}
        if not srcs:
            return [], {"n_source_clades": 0, "max_source_fraction": 1.0}
        for v in srcs.values():
            rng.shuffle(v)
        order = sorted(srcs, key=lambda s: (-len(srcs[s]), s))
        quota = -(-n // len(srcs))          # ceil
        picked, seen, taken = [], set(), defaultdict(int)
        # Round-robin so a large clade cannot fill the decoy before the small
        # ones are reached; the quota caps it even if everything else is empty.
        while len(picked) < n:
            progressed = False
            for s in order:
                if len(picked) >= n or taken[s] >= quota:
                    continue
                while srcs[s]:
                    k = srcs[s].pop()
                    if k in seen:
                        continue
                    seen.add(k)
                    picked.append(k)
                    taken[s] += 1
                    progressed = True
                    break
            if not progressed:
                break
        # Second pass: the quota is a budget, not a hard cap. Small source
        # clades run dry before N is reached, and refusing to fill would simply
        # delete the decoy -- shrinking the very denominator this null exists to
        # provide. Fill from whatever remains and let `max_source_fraction`
        # report honestly how concentrated the result ended up.
        if len(picked) < n:
            rest = [k for s in order for k in srcs[s] if k not in seen]
            rng.shuffle(rest)
            owner = {k: s for s in order for k in by_source[fam][s]}
            for k in rest:
                if len(picked) >= n:
                    break
                seen.add(k)
                picked.append(k)
                taken[owner.get(k, "?")] += 1
        frac = max(taken.values()) / len(picked) if picked else 1.0
        return picked, {"n_source_clades": sum(1 for s in taken if taken[s]),
                        "max_source_fraction": round(frac, 4),
                        "per_source": dict(sorted(taken.items()))}

    manifest = []
    skipped: dict[str, list[str]] = defaultdict(list)
    n_written = 0
    for fam, stem, n in sizes:
        avail = sum(len(v) for s, v in by_source[fam].items()
                    if args.allow_self_source or s != stem)
        if avail < n:
            print(f"  {stem}: {avail} windows outside the mirrored clade "
                  f"< clade size {n}, no null possible", file=sys.stderr)
            skipped[fam].append(stem)
            continue
        for rep in range(args.replicates):
            pick, diag = draw(fam, stem, n)
            if len(pick) < n:
                print(f"  {stem} rep{rep}: only {len(pick)}/{n} drawn, skipping",
                      file=sys.stderr)
                continue
            degenerate = (diag["n_source_clades"] < args.min_source_clades
                          or diag["max_source_fraction"] > args.max_source_fraction)
            name = f"DECOY_{fam}__rep{rep}__{stem.split('__', 1)[1]}"[:180]
            with open(args.out / f"{name}.fasta", "w") as fh:
                for k in pick:
                    h, s = pool[fam][k]
                    fh.write(f">{h}\n{s}\n")
            manifest.append({"decoy": name, "family": fam, "mirrors": stem,
                             "replicate": rep, "n_windows": n,
                             "null_strength": "degenerate" if degenerate else "strong",
                             **diag})
            n_written += 1

    # Per family: a null nobody can meet the quota for must not have its FDR
    # quoted as if it had been measured.
    per_family = {}
    for fam in sorted(pool):
        ds = [m for m in manifest if m["family"] == fam]
        if not ds:
            # A family whose clades admit no null at all is not "clean"; it is
            # unmeasured, and must say so rather than be absent from the table.
            per_family[fam] = {
                "decoys": 0, "degenerate": 0,
                "source_clades_available": len(by_source[fam]),
                "clades_with_no_null": sorted(skipped.get(fam, [])),
                "null_strength": "none"}
            print(f"  {fam}: no null possible "
                  f"({len(by_source[fam])} clade(s) in family)")
            continue
        bad = sum(1 for m in ds if m["null_strength"] == "degenerate")
        per_family[fam] = {
            "decoys": len(ds), "degenerate": bad,
            "clades_with_no_null": sorted(skipped.get(fam, [])),
            "source_clades_available": len(by_source[fam]),
            "median_source_clades": sorted(m["n_source_clades"] for m in ds)[len(ds) // 2],
            "max_source_fraction_worst": max(m["max_source_fraction"] for m in ds),
            "null_strength": "degenerate" if bad > len(ds) / 2 else "strong",
        }
        print(f"  {fam}: {len(ds)} decoys, {bad} degenerate, "
              f"null={per_family[fam]['null_strength']}")

    (args.out / "decoy_manifest.json").write_text(json.dumps(
        {"seed": args.seed, "replicates": args.replicates,
         "min_source_clades": args.min_source_clades,
         "max_source_fraction": args.max_source_fraction,
         "allow_self_source": args.allow_self_source,
         "real_clades": len(sizes), "decoy_clades": n_written,
         "per_family": per_family, "decoys": manifest}, indent=2))
    print(f"\n{n_written} decoy clades mirroring {len(sizes)} real clades -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
