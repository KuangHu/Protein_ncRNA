#!/usr/bin/env python
"""Locus sanity: which marker protein families sit within +/-N kb of each anchor.

Two uses, both from the census results:

  * **Cas3 verification.** PF18019 is an HD domain and the HD clan is
    promiscuous, so a Cas3 call is provisional until cascade genes (cas5/6/7/8,
    cas1/2) are seen beside it. An HD-only singleton is downweighted.
  * **RT subtyping.** Generic RT pools retrons, Abi-like defence RTs, DGRs,
    group II remnants and host enzymes. A neighbouring HNH/TOPRIM/HEPN effector
    distinguishes a defence RT from a bare one — without claiming 'retron',
    which cannot be called from RT homology alone.

Only the genomes carrying the selected anchors are re-processed, so this costs a
fraction of a full census.

    neighborhood_scan.py --census <corrected> --families Cas3 RT --out out/nbr
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import sys
import traceback
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_ncrna.compute import require_slurm
from protein_ncrna.anchors import build_anchor_hmm, find_anchors, validate_config
from protein_ncrna.orfs import EmptyGenome, call_orfs

_G: dict = {}


def _init(cfg, hmm, model_group, scratch, anchors_by_genome, flank):
    _G.update(cfg=cfg, hmm=hmm, model_group=model_group, scratch=scratch,
              anchors_by_genome=anchors_by_genome, flank=flank)


def _one(task: tuple[str, str]):
    """Return (genome, [(anchor_id, group, model, distance_bp)], error)."""
    acc, path = task
    try:
        orfs = call_orfs(path, workdir=_G["scratch"])
        hits = find_anchors(orfs, _G["hmm"], _G["model_group"], _G["cfg"],
                            workdir=_G["scratch"])
        if not hits:
            return acc, [], None
        by_contig = defaultdict(list)
        for h in hits:
            by_contig[h.orf.contig].append(h)

        out = []
        flank = _G["flank"]
        for aid, contig, a_start, a_end in _G["anchors_by_genome"].get(acc, ()):
            for h in by_contig.get(contig, ()):
                # Skip the anchor ORF itself.
                if h.orf.start == a_start and h.orf.end == a_end:
                    continue
                if h.orf.end < a_start - flank or h.orf.start > a_end + flank:
                    continue
                dist = (a_start - h.orf.end if h.orf.end < a_start
                        else h.orf.start - a_end if h.orf.start > a_end else 0)
                out.append((aid, h.family, h.model, dist))
        return acc, out, None
    except EmptyGenome as e:
        return acc, [], f"empty-stub: {e}"
    except Exception:
        return acc, [], traceback.format_exc(limit=3)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", required=True, type=Path)
    ap.add_argument("--genomes", type=Path,
                    default=Path("/global/scratch/users/kh36969/protein_ncrna/corpus/"
                                 "real_genomes.tsv"))
    ap.add_argument("--config", type=Path,
                    default=Path(__file__).resolve().parents[1] / "configs/blind/neighborhood.json")
    ap.add_argument("--families", nargs="+", required=True)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--scratch", default=None)
    ap.add_argument("--max-genomes", type=int, default=0,
                    help="cap the genomes scanned; a subsample answers the "
                         "neighbourhood question without the full corpus")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    args = ap.parse_args()
    require_slurm("neighborhood scan")
    args.out.mkdir(parents=True, exist_ok=True)

    # The neighbourhood config reuses the anchor machinery, so its "families"
    # are marker groups and `priority` is irrelevant.
    cfg = load_config_as_anchor_cfg(args.config)
    hmm = args.out / "neighborhood.hmm"
    model_group = build_anchor_hmm(cfg, hmm)
    flank = cfg["window"]["flank_bp"]
    print(f"[nbr] {len(model_group)} models in {len(cfg['families'])} groups, "
          f"+/-{flank} bp", file=sys.stderr)

    wanted = set(args.families)
    anchors_by_genome: dict[str, list] = defaultdict(list)
    with open(args.census / "anchors.tsv") as fh:
        cols = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[cols["family"]] not in wanted:
                continue
            aid = f[cols["anchor_id"]]
            # Recover anchor coordinates from the id: genome|contig:start-end:strand
            coord = aid.rsplit(":", 2)
            start, end = coord[-2].split("-")
            anchors_by_genome[f[cols["genome"]]].append(
                (aid, f[cols["contig"]], int(start), int(end)))
    print(f"[nbr] {sum(len(v) for v in anchors_by_genome.values())} anchors in "
          f"{len(anchors_by_genome)} genomes", file=sys.stderr)

    path_of = {}
    with open(args.genomes) as fh:
        cols = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            path_of[f[cols["accession"]]] = f[cols["path"]]

    todo = [(a, path_of[a]) for a in sorted(anchors_by_genome) if a in path_of]
    if args.max_genomes:
        todo = todo[::max(1, len(todo) // args.max_genomes)][:args.max_genomes]
    if args.nshards > 1:
        todo = todo[args.shard::args.nshards]
    print(f"[nbr] scanning {len(todo)} genomes on {args.workers} workers", file=sys.stderr)

    n_rows = 0
    fails = 0
    ctx = mp.get_context("fork")
    with open(args.out / "neighborhood.tsv", "w") as out:
        out.write("anchor_id\tneighbor_group\tneighbor_model\tdistance_bp\n")
        with ctx.Pool(args.workers, initializer=_init,
                      initargs=(cfg, str(hmm), model_group, args.scratch,
                                dict(anchors_by_genome), flank)) as pool:
            for i, (acc, rows, err) in enumerate(
                    pool.imap_unordered(_one, todo, chunksize=1), 1):
                if err:
                    fails += 1
                for aid, grp, model, dist in rows:
                    out.write(f"{aid}\t{grp}\t{model}\t{dist}\n")
                    n_rows += 1
                if i % 200 == 0 or i == len(todo):
                    # Flush with the progress line: a buffered handle means a
                    # job that hits its time limit loses every row it found.
                    out.flush()
                    print(f"[nbr] {i}/{len(todo)} genomes, {n_rows} neighbour hits",
                          file=sys.stderr, flush=True)
    print(f"[nbr] {n_rows} rows, {fails} failures -> {args.out}/neighborhood.tsv",
          file=sys.stderr)
    return 0


def load_config_as_anchor_cfg(path: Path) -> dict:
    """Adapt the neighbourhood config to the shape find_anchors expects.

    Marker groups become "families" so the same HMM search code runs; priority
    is uniform because a neighbour hitting two groups is genuinely ambiguous and
    both are worth recording.
    """
    import json as _json
    with open(path) as fh:
        raw = _json.load(fh)
    cfg = dict(raw)
    cfg["families"] = {g: {**spec, "priority": 10} for g, spec in raw["groups"].items()}
    cfg.pop("groups", None)
    cfg["window"] = {"flank_bp": raw.get("flank_bp", 10000),
                     "orient_to_anchor_strand": False, "min_contig_bp": 0}
    cfg["thresholds"] = {"use_cut_ga": True, "max_domain_evalue": 1e-5, "min_aa_len": 50}
    return validate_config(cfg)


if __name__ == "__main__":
    raise SystemExit(main())
