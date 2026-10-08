#!/usr/bin/env python
"""Arm 2: find more copies of a query protein across the genome corpus.

    expand_by_homology.py --genomes out/real_genomes.tsv --queries bag_reps.faa \\
                          --out <dir> --shard 0 --nshards 64

The insertion catalogue only extracts a protein when the same locus was observed
both with and without the element -- [genome][insertion][genome] against
[genome][genome]. That is a strong mobility proof and a severe sampling filter:
a copy is invisible unless a close relative happens to lack it. So a 6,964-site
bag and a 9-site bag may be equally common in nature and differ only in how often
an empty allele was sequenced.

This recovers the rest by homology. Each bag's representative proteins are
searched against every ORF in the corpus, and any ORF that matches becomes an
anchor in its own right, with no requirement that an empty allele exist. Those
anchors then carry the full +/-5 kb genomic window, where the bag carries only
the element -- a median of 249 bp upstream and 23 bp downstream, which truncates
any RNA extending past the insertion boundary.

Blind. The query proteins are MMseqs clusters of ab-initio Prodigal ORFs; the
search is sequence-only; no HMM, no IS library, no RNA model. A hit is a hit,
and what it is called is a question for after the predictions are frozen.

Structure mirrors collect_anchors.py: Prodigal per genome, shard over the corpus
so one failure costs one slice, and emit the same `anchors.tsv` schema so
extract_windows.py consumes the output unchanged. `family` carries the bag id of
the best-matching query, so each bag becomes its own clade downstream.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import subprocess
import sys
import tempfile
import traceback
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_ncrna.compute import require_slurm
from protein_ncrna.anchors import Anchor
from protein_ncrna.orfs import EmptyGenome, call_orfs
from protein_ncrna.seqio import read_fasta
from protein_ncrna.tools import which
from protein_ncrna.windows import flank_availability

ANCHOR_COLS = ["anchor_id", "genome", "species", "organism", "contig",
               "contig_len", "family", "model", "all_models", "all_families",
               "ambiguous", "score", "evalue", "anchor_strand", "anchor_len",
               "anchor_partial", "win_start", "win_end", "anchor_offset",
               "left_flank_bp", "right_flank_bp", "truncated_left",
               "truncated_right", "window_complete"]

WCFG = {"window": {"flank_bp": 5000, "orient_to_anchor_strand": True,
                   "min_contig_bp": 11000}}

_G: dict = {}


def _init(qdb, scratch, species, min_id, min_cov, evalue):
    _G.update(qdb=qdb, scratch=scratch, species=species, min_id=min_id,
              min_cov=min_cov, evalue=evalue)


def _search(acc: str, path: str) -> list[dict]:
    """ORFs of one genome that match any query protein."""
    try:
        orfs = call_orfs(path, workdir=_G["scratch"])
    except EmptyGenome:
        return []
    except Exception:
        traceback.print_exc()
        return []
    if not orfs:
        return []
    lens: dict[str, int] = {}
    for h, s in read_fasta(path):
        lens[h.split()[0]] = len(s)

    with tempfile.TemporaryDirectory(dir=_G["scratch"]) as td:
        td = Path(td)
        q = td / "orfs.faa"
        with open(q, "w") as fh:
            for i, o in enumerate(orfs):
                fh.write(f">{i}\n{o.aa}\n")
        res = td / "hits.m8"
        proc = subprocess.run(
            [which("mmseqs"), "easy-search", str(q), str(_G["qdb"]), str(res),
             str(td / "t"), "--threads", "1", "-v", "1",
             "-e", str(_G["evalue"]), "--max-seqs", "20", "-s", "5.7",
             "--format-output",
             "query,target,fident,alnlen,evalue,bits,qlen,tlen"],
            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        if proc.returncode != 0 or not res.exists():
            return []
        best: dict[int, tuple] = {}
        allq: dict[int, set] = defaultdict(set)
        for line in open(res):
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            i = int(f[0])
            fid, aln, ev, bits = float(f[2]), int(f[3]), float(f[4]), float(f[5])
            qlen, tlen = int(f[6]), int(f[7])
            # Coverage on the *query* protein, so a short ORF matching one domain
            # of a long bag protein is not recruited as a copy of it.
            if fid < _G["min_id"] or aln / max(1, tlen) < _G["min_cov"]:
                continue
            allq[i].add(f[1].split("|")[0])
            if i not in best or bits > best[i][1]:
                best[i] = (f[1], bits, ev, fid)

    rows = []
    for i, (tgt, bits, ev, fid) in sorted(best.items()):
        o = orfs[i]
        bag = tgt.split("|")[0]
        geom = flank_availability(Anchor(o, bag, tgt, ev, bits),
                                  lens.get(o.contig, 0), WCFG)
        if geom is None:
            continue
        fams = sorted(allq[i])
        rows.append({
            "anchor_id": f"{acc}|{o.contig}:{o.start}-{o.end}:"
                         f"{'+' if o.strand > 0 else '-'}",
            "genome": acc, "species": _G["species"].get(acc, "unknown"),
            "organism": _G["species"].get(acc, "unknown"),
            "contig": o.contig, "contig_len": lens.get(o.contig, 0),
            "family": bag, "model": tgt,
            "all_models": ",".join(fams), "all_families": ",".join(fams),
            "ambiguous": len(fams) > 1,
            "score": round(bits, 1), "evalue": ev,
            "anchor_strand": o.strand, "anchor_len": o.end - o.start + 1,
            "anchor_partial": o.partial, "aa": o.aa,
            **geom,
        })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes", required=True, type=Path)
    ap.add_argument("--queries", required=True, type=Path,
                    help="FASTA of bag representative proteins, ids '<bag>|<n>'")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--scratch", default=None)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--min-identity", type=float, default=0.30)
    ap.add_argument("--min-coverage", type=float, default=0.70)
    ap.add_argument("--evalue", type=float, default=1e-10)
    args = ap.parse_args()
    require_slurm("homology expansion")
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "proteins").mkdir(exist_ok=True)
    scratch = args.scratch or tempfile.mkdtemp()

    rows, species = [], {}
    with open(args.genomes) as fh:
        cols = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            acc, p = f[cols["accession"]], f[cols["path"]]
            rows.append((acc, p))
            if "species" in cols and len(f) > cols["species"]:
                species[acc] = f[cols["species"]]
    rows.sort()
    if args.limit:
        rows = rows[:args.limit]
    mine = rows[args.shard::args.nshards]
    print(f"shard {args.shard}/{args.nshards} of {len(rows)} genomes "
          f"-> {len(mine)} to process", flush=True)

    nq = sum(1 for _ in read_fasta(str(args.queries)))
    print(f"{nq} query proteins", flush=True)

    out_rows: list[dict] = []
    with mp.Pool(args.workers, initializer=_init,
                 initargs=(args.queries, scratch, species, args.min_identity,
                           args.min_coverage, args.evalue)) as pool:
        for n, got in enumerate(pool.starmap(_search, mine, chunksize=4), 1):
            out_rows.extend(got)
            if n % 500 == 0:
                print(f"  {n}/{len(mine)} genomes, {len(out_rows)} anchors",
                      flush=True)

    tsv = args.out / f"anchors_shard{args.shard}.tsv"
    with open(tsv, "w") as fh:
        fh.write("\t".join(ANCHOR_COLS) + "\n")
        for r in out_rows:
            fh.write("\t".join(str(r[c]) for c in ANCHOR_COLS) + "\n")
    byfam: dict[str, list] = defaultdict(list)
    for r in out_rows:
        byfam[r["family"]].append(r)
    for fam, rs in byfam.items():
        with open(args.out / "proteins" / f"{fam}.shard{args.shard}.faa", "w") as fh:
            for r in rs:
                fh.write(f">{r['anchor_id']}\n{r['aa']}\n")

    summary = {
        "shard": args.shard, "nshards": args.nshards,
        "genomes_processed": len(mine), "anchors": len(out_rows),
        "complete_windows": sum(1 for r in out_rows if r["window_complete"]),
        "bags_hit": len(byfam),
        "min_identity": args.min_identity, "min_coverage": args.min_coverage,
        "evalue": args.evalue,
    }
    (args.out / f"summary_shard{args.shard}.json").write_text(
        json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"-> {tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
