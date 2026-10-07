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
        kept: dict[str, int] = defaultdict(int)
        for acc in list(by_genome):
            keep = []
            for w in by_genome[acc]:
                if any(kept[g] < args.max_per_clade for g in w["groups"]):
                    for g in w["groups"]:
                        kept[g] += 1
                    keep.append(w)
            by_genome[acc] = keep

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
