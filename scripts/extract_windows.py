#!/usr/bin/env python
"""Carve the oriented +/-5 kb windows for a chosen set of anchors.

The census deliberately stored no sequence, so windows are cut on demand. This
needs no ORF calling: `anchors.tsv` already records `win_start`, `win_end` and
`anchor_strand`, so a window is a slice of the contig plus a reverse complement,
reproducing `windows.extract_window` exactly. That makes it IO-bound and far
cheaper than the census.

    # per selected clade, for step 4
    extract_windows.py --census <corrected> --clades <clades>/selected_clades.tsv \\
                       --members <clades>/cluster_members.tsv --out out/windows

    # or a whole family, for the retron benchmark
    extract_windows.py --census <corrected> --family RT --out out/rt_windows

Only complete windows are written by default: a truncated window cannot support
a conserved-block comparison and would silently weaken any alignment.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_ncrna.compute import require_slurm
from protein_ncrna.diversity import (farthest_point, hash_sample,
                                     stratified_sample)
from protein_ncrna.seqio import read_fasta, revcomp

_G: dict = {}


def _init(by_genome, keep_incomplete):
    _G.update(by_genome=by_genome, keep_incomplete=keep_incomplete)


def _one(task):
    acc, path = task
    want = _G["by_genome"][acc]
    try:
        need = {w["contig"] for w in want}
        contigs = {}
        for h, s in read_fasta(path):
            cid = h.split()[0]
            if cid in need:
                contigs[cid] = s
        out = []
        for w in want:
            seq = contigs.get(w["contig"])
            if seq is None:
                continue
            sub = seq[w["win_start"] - 1:w["win_end"]]
            if len(sub) != w["win_end"] - w["win_start"] + 1:
                continue
            if w["strand"] < 0:
                sub = revcomp(sub)
            out.append({**w, "seq": sub})
        return acc, out, None
    except Exception as e:
        return acc, [], f"{type(e).__name__}: {e}"


def load_family_proteins(census: Path, families: set[str]) -> dict[str, str]:
    """Anchor protein sequences, keyed by anchor id, for the given families."""
    seqs: dict[str, str] = {}
    for fam in sorted(families):
        faa = census / "proteins" / f"{fam}.faa"
        if not faa.exists():
            print(f"  no {faa}; cap falls back to a hash sample for {fam}",
                  file=sys.stderr)
            continue
        for h, s in read_fasta(str(faa)):
            seqs.setdefault(h.split()[0], s)
    return seqs


def select_capped(by_genome: dict, args) -> list[tuple]:
    """Apply --max-per-clade per clade and prune `by_genome` in place.

    Returns the manifest rows (clade_id, anchor_id, selected, reason) for every
    candidate anchor, kept or dropped, so the cap is auditable rather than
    implicit in which windows happen to exist.

    One anchor can belong to several clades (the thresholds are nested), so an
    anchor is written if *any* of its clades selected it; `reason` records the
    clade that this row is about.
    """
    per_clade: dict[str, list[str]] = defaultdict(list)
    win: dict[str, dict] = {}
    fams: set[str] = set()
    for ws in by_genome.values():
        for w in ws:
            win[w["anchor_id"]] = w
            fams.add(w["family"])
            for g in w["groups"]:
                per_clade[g].append(w["anchor_id"])

    seqs = (load_family_proteins(args.census, fams)
            if args.cap_mode in ("stratified", "diverse") else {})
    if args.cap_mode in ("stratified", "diverse") and not seqs:
        # Without proteins both modes silently degrade to a hash sample, which
        # is order-independent but not diversity-aware -- exactly the thing the
        # cap exists to provide. Say so rather than letting the manifest's
        # reason column be the only trace.
        print(f"WARNING: no anchor proteins under {args.census / 'proteins'}; "
              f"--cap-mode {args.cap_mode} degrades to a hash sample",
              file=sys.stderr)
    # Prefer complete windows, then longer ones: a truncated window cannot
    # support a block comparison, so it should lose a tie even if it is the
    # more diverse choice.
    prefer = {a: (bool(w["window_complete"]), w["win_end"] - w["win_start"])
              for a, w in win.items()}

    rows, keep_ids = [], set()
    for clade in sorted(per_clade):
        ids = sorted(set(per_clade[clade]))
        if len(ids) <= args.max_per_clade:
            sel, why = set(ids), "under_cap"
        elif args.cap_mode == "first":
            sel, why = set(ids[:args.max_per_clade]), "input_order"
        elif args.cap_mode == "hash":
            sel = set(hash_sample(ids, args.max_per_clade, args.seed))
            why = "hash_sample"
        else:
            have = sum(1 for i in ids if i in seqs)
            if have < len(ids) * 0.5:
                sel = set(hash_sample(ids, args.max_per_clade, args.seed))
                why = f"hash_sample_fallback_{have}_of_{len(ids)}_proteins"
            elif args.cap_mode == "diverse":
                sel = set(farthest_point(ids, seqs, args.max_per_clade,
                                         seed=args.seed, prefer=prefer))
                why = "diverse_kmer"
            else:
                # Stratified, not farthest-point. Taking the 500 mutually most
                # dissimilar members of a 7,877-member clade does not sample its
                # breadth -- it selects its 500 outliers, which is measurably
                # worse than the arbitrary slice it replaced: on this corpus it
                # cut GroupII_RT from 4 real blocks to 2 while decoys rose to 13.
                # Seeding farthest-point and then keeping each seed's *nearest*
                # representative covers the same breadth with typical members.
                sel = set(stratified_sample(ids, seqs, args.max_per_clade,
                                            seed=args.seed, prefer=prefer))
                why = "stratified_kmer"
        keep_ids |= sel
        for i in ids:
            rows.append((clade, i, str(i in sel), why))
        if len(ids) > args.max_per_clade:
            print(f"  cap {clade[:58]}: {len(ids)} -> {len(sel)} ({why})",
                  flush=True)

    for acc in list(by_genome):
        by_genome[acc] = [w for w in by_genome[acc] if w["anchor_id"] in keep_ids]
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", required=True, type=Path)
    ap.add_argument("--genomes", type=Path,
                    default=Path("/global/scratch/users/kh36969/protein_ncrna/corpus/"
                                 "real_genomes.tsv"))
    ap.add_argument("--clades", type=Path, help="selected_clades.tsv")
    ap.add_argument("--members", type=Path, help="cluster_members.tsv")
    ap.add_argument("--family", help="extract a whole family instead of clades")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--keep-incomplete", action="store_true")
    ap.add_argument("--proteins", action="store_true",
                    help="also write one anchor-protein FASTA per clade, from "
                         "the census per-family FASTAs (the step-3 tree input)")
    ap.add_argument("--max-per-clade", type=int, default=0,
                    help="cap anchors per clade; a clade of 10k near-identical "
                         "copies does not need every one aligned")
    ap.add_argument("--cap-mode",
                    choices=["stratified", "diverse", "hash", "first"],
                    default="stratified",
                    help="how the cap chooses which anchors to keep: spread "
                         "across anchor-protein diversity (default), a "
                         "deterministic order-independent sample, or input "
                         "order (the old behaviour; reproduces earlier runs)")
    ap.add_argument("--seed", type=int, default=1,
                    help="tie-break seed for --cap-mode diverse/hash")
    args = ap.parse_args()
    require_slurm("window extraction")
    args.out.mkdir(parents=True, exist_ok=True)

    # anchor -> which clade(s) it belongs to
    group_of: dict[str, list[str]] = defaultdict(list)
    if args.clades and args.members:
        sel = set()
        with open(args.clades) as fh:
            c = {k: i for i, k in enumerate(next(fh).rstrip("\n").split("\t"))}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                sel.add((f[c["family"]], f[c["threshold"]], f[c["cluster_id"]]))
        with open(args.members) as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if (f[0], f[1], f[2]) in sel:
                    group_of[f[3]].append(f"{f[0]}__thr{f[1]}__{_safe(f[2])}")
        print(f"{len(sel)} selected clades, {len(group_of)} member anchors")
    elif not args.family:
        raise SystemExit("give either --clades with --members, or --family")

    by_genome: dict[str, list] = defaultdict(list)
    n_total = n_incomplete = 0
    with open(args.census / "anchors.tsv") as fh:
        c = {k: i for i, k in enumerate(next(fh).rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            aid = f[c["anchor_id"]]
            if args.family:
                if f[c["family"]] != args.family:
                    continue
                groups = [args.family]
            else:
                groups = group_of.get(aid)
                if not groups:
                    continue
            n_total += 1
            if f[c["window_complete"]] != "True":
                n_incomplete += 1
                if not args.keep_incomplete:
                    continue
            by_genome[f[c["genome"]]].append({
                "anchor_id": aid, "groups": groups,
                "contig": f[c["contig"]],
                "win_start": int(f[c["win_start"]]), "win_end": int(f[c["win_end"]]),
                "strand": int(f[c["anchor_strand"]]),
                "anchor_offset": int(f[c["anchor_offset"]]),
                "anchor_len": int(f[c["anchor_len"]]),
                "family": f[c["family"]], "species": f[c["species"]],
                "window_complete": f[c["window_complete"]] == "True",
            })

    if args.max_per_clade:
        # The cap decides which half of a 9,600-member clade step 4 ever sees,
        # so it must not be "the first 500 in census order". That order tracks
        # manifest and accession order, which correlates with submission batch
        # and therefore with exactly the subclade structure the cap should be
        # spreading across -- and when R-scape later reports a block as
        # underpowered there is no way to tell a biologically tight clade from
        # a narrowly sampled one.
        cap_rows = select_capped(by_genome, args)
        with open(args.out / "cap_manifest.tsv", "w") as fh:
            fh.write("clade_id\tanchor_id\tselected\treason\n")
            for r in cap_rows:
                fh.write("\t".join(r) + "\n")

    # Anchor proteins per clade. Written from the same post-cap membership as
    # the windows, so the protein set and the window set are the same anchors —
    # a clade tree built on a different subset than the alignment is misleading.
    if args.proteins:
        pdir = args.out / "proteins"
        pdir.mkdir(exist_ok=True)
        want: dict[str, set[str]] = defaultdict(set)
        fam_of: dict[str, str] = {}
        for ws in by_genome.values():
            for w in ws:
                for g in w["groups"]:
                    want[g].add(w["anchor_id"])
                    fam_of[g] = w["family"]
        n_p = 0
        for fam in sorted(set(fam_of.values())):
            faa = args.census / "proteins" / f"{fam}.faa"
            if not faa.exists():
                print(f"  no {faa}, skipping proteins for {fam}", file=sys.stderr)
                continue
            groups = [g for g, f in fam_of.items() if f == fam]
            seqs = {h.split()[0]: (h, s) for h, s in read_fasta(str(faa))}
            for g in groups:
                with open(pdir / f"{g}.faa", "w") as fh:
                    for aid in sorted(want[g]):
                        if aid in seqs:
                            h, s = seqs[aid]
                            fh.write(f">{h}\n{s}\n")
                            n_p += 1
        print(f"{n_p} anchor proteins in {len(want)} clade FASTAs -> {pdir}")

    n_want = sum(len(v) for v in by_genome.values())
    print(f"{n_want} windows to carve from {len(by_genome)} genomes "
          f"({n_incomplete}/{n_total} anchors dropped as incomplete)")

    path_of = {}
    with open(args.genomes) as fh:
        c = {k: i for i, k in enumerate(next(fh).rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            path_of[f[c["accession"]]] = f[c["path"]]
    todo = [(a, path_of[a]) for a in sorted(by_genome) if a in path_of]

    handles: dict[str, object] = {}
    fastas: dict[str, object] = {}
    n_written = fails = 0
    ctx = mp.get_context("fork")
    with ctx.Pool(args.workers, initializer=_init,
                  initargs=(dict(by_genome), args.keep_incomplete)) as pool:
        for i, (acc, recs, err) in enumerate(pool.imap_unordered(_one, todo, chunksize=4), 1):
            if err:
                fails += 1
            for r in recs:
                for g in r["groups"]:
                    if g not in handles:
                        handles[g] = open(args.out / f"{g}.jsonl", "w")
                        fastas[g] = open(args.out / f"{g}.fasta", "w")
                    rec = {k: v for k, v in r.items() if k != "groups"}
                    handles[g].write(json.dumps(rec) + "\n")
                    fastas[g].write(
                        f">{r['anchor_id']} family={r['family']} species={r['species']} "
                        f"anchor_offset={r['anchor_offset']} anchor_len={r['anchor_len']}\n"
                        f"{r['seq']}\n")
                    n_written += 1
            if i % 500 == 0 or i == len(todo):
                for h in handles.values():
                    h.flush()
                for h in fastas.values():
                    h.flush()
                print(f"  {i}/{len(todo)} genomes, {n_written} windows", flush=True)
    for h in list(handles.values()) + list(fastas.values()):
        h.close()

    print(f"{n_written} windows in {len(handles)} groups, {fails} genome failures "
          f"-> {args.out}")
    return 0


def _safe(s: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in s)[:80]


if __name__ == "__main__":
    raise SystemExit(main())
