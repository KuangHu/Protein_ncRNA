#!/usr/bin/env python
"""Decoy sets for Arm 2, where the frozen within-family permutation cannot run.

    arm2_nulls.py --windows <pilot_windows/BAG> --mode genomic_background \\
                  --genomes out/real_genomes.tsv --out <decoys> --replicates 3

The frozen null substitutes windows drawn from *other clades of the same
family*. Arm 2 has no such family for the bags that matter: 124 of 399 are
singletons, and the four most diverse bags -- including CDS05629 with 28,896
distinct proteins -- have no sibling at all. A null with no source clades still
prints an FDR of 0.0, which is the failure that invalidated Cas3 and emptied
every Arm 1 FDR cell, so the answer is a different null rather than a missing
one.

Two are built here, and they ask different questions.

`genomic_background` draws random windows of the same length from the same
genomes the real bag's members came from. It asks: would *any* neighbourhood in
these genomes produce a block like this? It is available for every bag, and it
is an easier null than sibling permutation -- unrelated windows share less
sequence, so real-vs-decoy separation is overstated rather than matched. It
bounds the answer; it does not measure it, and must be labelled that way.

`within_bag_shuffle` keeps the bag's own sequences and moves the anchor. Same
windows, same composition, same relatedness -- only the anchor frame is
randomised, kept at least `--min-anchor-sep` from the true anchor. It asks: is
the recurrence *anchored*, or would any frame on these sequences recover it?
That is the strictest of the three, because every source of recurrence except
the anchor's position survives into the decoy.

Blind. Only coordinates, orientation and genomic sequence are touched; no RNA
model, length prior or curated element enters here.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import json
import random
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_ncrna.compute import require_slurm
from protein_ncrna.seqio import read_fasta, revcomp

# A pseudo-anchor placed next to the real one would sit in the same noncoding
# run and recover the same block, which would make the shuffle null look
# identical to the real set for reasons that have nothing to do with chance.
MIN_ANCHOR_SEP = 2000


def load_windows(d: Path) -> tuple[str, list[dict]]:
    """The one clade in a pilot window dir, and its records."""
    jl = sorted(d.glob("*.jsonl"))
    if len(jl) != 1:
        raise SystemExit(f"expected exactly one clade jsonl in {d}, found {len(jl)}")
    clade = jl[0].stem
    recs = [json.loads(line) for line in open(jl[0])]
    if not recs:
        raise SystemExit(f"{jl[0]} is empty")
    return clade, recs


def write_clade(out: Path, clade: str, recs: list[dict]) -> None:
    out.mkdir(parents=True, exist_ok=True)
    with open(out / f"{clade}.fasta", "w") as fa, open(out / f"{clade}.jsonl", "w") as jl:
        for r in recs:
            jl.write(json.dumps(r) + "\n")
            fa.write(f">{r['anchor_id']} family={r['family']} "
                     f"species={r.get('species', '')} "
                     f"anchor_offset={r['anchor_offset']} "
                     f"anchor_len={r['anchor_len']}\n{r['seq']}\n")


def write_shuffled_proteins(pin: Path, pout: Path, clade: str,
                            recs: list[dict], rep: int) -> int:
    """Carry each shuffled window's anchor protein across under its new id.

    discover_blocks picks references by anchor-protein 5-mer diversity and falls
    back, per id, to the window's first 3 kb of DNA when the protein is missing.
    The shuffle decoys rename every anchor to `<id>__sh<rep>`, so without this
    they would miss that lookup entirely and pick references in DNA space while
    the real arm picks in protein space. That is a second difference on top of
    the one the pilot is trying to measure, and it lands hardest on the arm that
    is supposed to be the cleanest: same windows, same proteins, only the anchor
    frame moved.
    """
    src = pin / f"{clade}.faa"
    if not src.exists():
        print(f"  no {src}; shuffled decoys will fall back to DNA references",
              file=sys.stderr)
        return 0
    seqs = {h.split()[0]: s for h, s in read_fasta(str(src))}
    pout.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(pout / f"DECOY_within_bag_shuffle__rep{rep}__{clade}.faa", "w") as fh:
        for r in recs:
            base = r["anchor_id"].rsplit("__sh", 1)[0]
            s = seqs.get(base)
            if s:
                fh.write(f">{r['anchor_id']}\n{s}\n")
                n += 1
    return n


def shuffle_anchor(recs: list[dict], rng: random.Random, rep: int,
                   min_sep: int) -> list[dict]:
    """Same sequences, anchor frame moved and orientation randomised."""
    out = []
    for r in recs:
        seq = r["seq"]
        alen = int(r["anchor_len"])
        if len(seq) < alen + 2 * min_sep:
            continue
        real = int(r["anchor_offset"])
        # Offsets at least min_sep from the real anchor, on either side.
        lo_hi = (0, real - min_sep - alen)
        hi_lo = (real + alen + min_sep, len(seq) - alen)
        spans = [s for s in (lo_hi, hi_lo) if s[1] > s[0]]
        if not spans:
            continue
        a, b = rng.choice(spans)
        off = rng.randint(a, b)
        # Both offsets must be flipped together. Reporting the pseudo-anchor in
        # the reverse-complemented frame and the real one in the original frame
        # would silently corrupt `|pseudo - real|` for half the records, which
        # is the one audit this field exists to support.
        real_out = real
        if rng.random() < 0.5:
            seq = revcomp(seq)
            off = len(seq) - off - alen
            real_out = len(seq) - real - alen
        out.append({**r, "anchor_id": f"{r['anchor_id']}__sh{rep}",
                    "seq": seq, "anchor_offset": off, "anchor_len": alen,
                    "pseudo_anchor": True, "real_anchor_offset": real_out,
                    "real_anchor_offset_input_frame": real})
    return out


def genomic_background(recs: list[dict], genomes: dict[str, str],
                       rng: random.Random, replicates: int) -> list[list[dict]]:
    """Random same-length windows from the same genomes, random orientation.

    Every replicate is produced in the same pass over the genomes. Reading each
    gzipped genome once per replicate instead would multiply the only expensive
    step in this script by `replicates` for no benefit -- the draws are
    independent either way.
    """
    by_genome: dict[str, list[dict]] = defaultdict(list)
    for r in recs:
        g = r.get("genome") or r["anchor_id"].split("|")[0]
        by_genome[g].append(r)

    out: list[list[dict]] = [[] for _ in range(replicates)]
    n_missing = n_unreadable = 0
    for g, group in sorted(by_genome.items()):
        path = genomes.get(g)
        if path is None:
            n_missing += 1
            continue
        try:
            contigs = [(h.split()[0], s) for h, s in read_fasta(path) if s]
        except (OSError, EOFError, gzip.BadGzipFile):
            n_unreadable += 1
            continue
        if not contigs:
            n_unreadable += 1
            continue
        contigs.sort(key=lambda c: len(c[1]))
        lens = [len(s) for _, s in contigs]
        # Window lengths are near-constant within a bag, so the set of contigs
        # long enough to hold one is worth resolving once per length rather than
        # once per record.
        usable_from: dict[int, int] = {}
        for rep in range(replicates):
            for i, r in enumerate(group):
                wlen = len(r["seq"])
                alen = int(r["anchor_len"])
                lo = usable_from.get(wlen)
                if lo is None:
                    lo = bisect.bisect_left(lens, wlen)
                    usable_from[wlen] = lo
                if lo >= len(contigs):
                    continue
                name, cseq = contigs[rng.randrange(lo, len(contigs))]
                start = rng.randint(0, len(cseq) - wlen)
                seq = cseq[start:start + wlen].upper()
                off = rng.randint(0, max(0, wlen - alen))
                if rng.random() < 0.5:
                    seq = revcomp(seq)
                    off = len(seq) - off - alen
                out[rep].append({
                    "anchor_id": f"BG{rep}_{g}_{i}", "contig": name,
                    "genome": g, "species": r.get("species", ""),
                    "win_start": start, "win_end": start + wlen,
                    "strand": 1, "anchor_offset": off, "anchor_len": alen,
                    "family": r["family"], "window_complete": True,
                    "pseudo_anchor": True, "seq": seq})
    if n_missing or n_unreadable:
        print(f"  {n_missing} genomes not in the table, "
              f"{n_unreadable} unreadable", file=sys.stderr)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", required=True, type=Path,
                    help="one bag's pilot window dir")
    ap.add_argument("--mode", required=True,
                    choices=["genomic_background", "within_bag_shuffle"])
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--genomes", type=Path,
                    help="accession->path tsv; required for genomic_background")
    ap.add_argument("--replicates", type=int, default=3)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--min-anchor-sep", type=int, default=MIN_ANCHOR_SEP)
    ap.add_argument("--proteins-in", type=Path,
                    help="anchor-protein dir for the real clade "
                         "(default <windows>/proteins); within_bag_shuffle "
                         "re-emits these under the shuffled anchor ids")
    ap.add_argument("--proteins-out", type=Path,
                    help="where to write the shuffled anchor proteins "
                         "(default <out>/proteins)")
    args = ap.parse_args()
    require_slurm("arm2 null generation")

    clade, recs = load_windows(args.windows)
    rng = random.Random(f"{args.seed}:{args.mode}:{clade}")
    print(f"{clade}: {len(recs)} real windows, mode={args.mode}")

    genomes: dict[str, str] = {}
    if args.mode == "genomic_background":
        if not args.genomes:
            raise SystemExit("--genomes is required for genomic_background")
        with open(args.genomes) as fh:
            c = {k: i for i, k in enumerate(next(fh).rstrip("\n").split("\t"))}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                genomes[f[c["accession"]]] = f[c["path"]]

    if args.mode == "within_bag_shuffle":
        reps = [shuffle_anchor(recs, rng, rep, args.min_anchor_sep)
                for rep in range(args.replicates)]
    else:
        reps = genomic_background(recs, genomes, rng, args.replicates)

    pin = args.proteins_in or (args.windows / "proteins")
    pout = args.proteins_out or (args.out / "proteins")

    made = []
    for rep, d in enumerate(reps):
        name = f"DECOY_{args.mode}__rep{rep}__{clade}"
        if len(d) < 4:
            print(f"  {name}: only {len(d)} windows, skipped", file=sys.stderr)
            continue
        write_clade(args.out, name, d)
        np = (write_shuffled_proteins(pin, pout, clade, d, rep)
              if args.mode == "within_bag_shuffle" else 0)
        made.append({"clade": name, "n_windows": len(d),
                     "n_anchor_proteins": np,
                     "retained_frac": round(len(d) / len(recs), 4)})
        print(f"  {name}: {len(d)} windows ({len(d)/len(recs):.1%} of real)"
              + (f", {np} anchor proteins" if args.mode == "within_bag_shuffle"
                 else ""))

    # A decoy set far smaller than the real one is not a matched null: block
    # recurrence is a fraction of members, so a short decoy clade clears the
    # recurrence floor less often for reasons of size alone.
    (args.out / f"null_manifest.{args.mode}.{clade}.json").write_text(json.dumps({
        "clade": clade, "mode": args.mode, "n_real": len(recs),
        "replicates": args.replicates, "seed": args.seed,
        "min_anchor_sep": args.min_anchor_sep if
        args.mode == "within_bag_shuffle" else None,
        "decoys": made,
        # Which space discover_blocks will pick this arm's references in. The
        # genomic background has no anchor protein by construction, so it is an
        # explicit, declared fallback rather than an accident -- and one more
        # reason its separation reads as a bound, not a match.
        "reference_space": ("sequence_fallback"
                            if args.mode == "genomic_background"
                            else "protein" if made and
                            all(m["n_anchor_proteins"] for m in made)
                            else "sequence_fallback"),
        "null_strength": ("matched" if made and
                          min(m["retained_frac"] for m in made) >= 0.9
                          else "short" if made else "none"),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
