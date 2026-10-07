#!/usr/bin/env python
"""Steps 1-2: genomes -> anchor proteins -> (optionally) oriented +/-5 kb windows.

Two modes:

    # census: anchors + flank feasibility + per-family protein FASTA, no sequence
    collect_anchors.py --genomes real_genomes.tsv --out out/census --windows none

    # discovery: additionally carve and store every window
    collect_anchors.py --genomes real_genomes.tsv --out out/disc --windows all

`--genomes` is a TSV with `accession` and `path` columns (`species` optional) —
use `out/real_genomes.tsv` from genomes_db_report.py. Writes:

    <out>/anchors.tsv                   one row per anchor, with flank feasibility
    <out>/proteins/<family>.faa         anchor protein sequences, for divergence
    <out>/summary_by_anchor_family.tsv  per-family counts
    <out>/summary.json
    <out>/windows.jsonl                 only when --windows all

Window sequence is ~10 kb per anchor; at census scale that is tens of GB for
data the census never reads, which is why `none` is the default.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
import traceback
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_ncrna.compute import require_slurm
from protein_ncrna.anchors import build_anchor_hmm, find_anchors, load_config
from protein_ncrna.orfs import EmptyGenome, call_orfs
from protein_ncrna.seqio import load_contigs, read_fasta
from protein_ncrna.windows import extract_window, flank_availability

_G: dict = {}


def _init(cfg, anchor_hmm, model_family, scratch, want_seq):
    _G.update(cfg=cfg, anchor_hmm=anchor_hmm, model_family=model_family,
              scratch=scratch, want_seq=want_seq)


def _contig_lengths(path: str) -> tuple[dict[str, int], str]:
    """Contig lengths without holding the sequences in memory, plus organism.

    40% of the real genomes in this mirror sit outside the species views and
    would otherwise be recorded as 'unknown'. The assembly's own FASTA header
    names the organism, and the file is being read anyway.
    """
    lengths, first = {}, ""
    for h, s in read_fasta(path):
        if not first:
            first = h
        lengths[h.split()[0]] = len(s)
    return lengths, organism_from_header(first)


def organism_from_header(header: str) -> str:
    """'CP000026.1 Salmonella enterica subsp. ...' -> 'salmonella_enterica'.

    Returned in the same slug form the species views use, so a genome labelled
    by a view and one labelled from its header collapse to a single species
    rather than being counted twice in the diversity summary.
    """
    parts = header.split()[1:]
    if len(parts) < 2:
        return "unknown"
    genus, species = parts[0].strip(","), parts[1].strip(",")
    if not genus[:1].isupper() or not species[:1].islower() or not species.isalpha():
        return "unknown"
    return f"{genus.lower()}_{species.lower()}"


def _one(task: tuple[str, str]) -> tuple[str, list[dict], list[dict], str | None]:
    """Return (accession, anchor rows, window records, error)."""
    accession, path = task
    cfg, want_seq = _G["cfg"], _G["want_seq"]
    try:
        orfs = call_orfs(path, workdir=_G["scratch"])
        anchors = find_anchors(orfs, _G["anchor_hmm"], _G["model_family"], cfg,
                               workdir=_G["scratch"])
        if not anchors:
            return accession, [], [], None

        if want_seq:
            # One pass: keep sequences, lengths, and the first full header.
            # load_contigs() drops the description, which is where the organism is.
            contigs, lengths, first = {}, {}, ""
            for h, s in read_fasta(path):
                if not first:
                    first = h
                cid = h.split()[0]
                contigs[cid] = s
                lengths[cid] = len(s)
            organism = organism_from_header(first)
        else:
            contigs = None
            lengths, organism = _contig_lengths(path)

        rows, windows = [], []
        for a in anchors:
            clen = lengths.get(a.orf.contig)
            if clen is None:
                continue
            geom = flank_availability(a, clen, cfg)
            if geom is None:
                continue
            rows.append({
                "anchor_id": f"{accession}|{a.orf.orf_id}",
                "genome": accession,
                "contig": a.orf.contig,
                "contig_len": clen,
                "family": a.family,
                "model": a.model,
                "all_models": ",".join(a.all_models),
                "all_families": ",".join(a.all_families),
                "ambiguous": a.is_ambiguous,
                "score": a.score,
                "evalue": a.evalue,
                "anchor_strand": a.orf.strand,
                "anchor_len": a.orf.end - a.orf.start + 1,
                "anchor_partial": a.orf.is_partial,
                "aa": a.orf.aa,
                **geom,
            })
            if want_seq:
                w = extract_window(a, contigs[a.orf.contig], accession, cfg)
                if w is not None:
                    windows.append(w.to_record())
        for r in rows:
            r["organism"] = organism
        return accession, rows, windows, None
    except EmptyGenome as e:
        return accession, [], [], f"empty-stub: {e}"
    except Exception:
        return accession, [], [], traceback.format_exc(limit=3)


def read_genomes(path: Path) -> tuple[list[tuple[str, str]], dict[str, str]]:
    rows, species = [], {}
    with open(path) as fh:
        first = fh.readline().rstrip("\n").split("\t")
        cols = {c: i for i, c in enumerate(first)}
        if "accession" not in cols:  # headerless
            fh.seek(0)
            cols = {"accession": 0, "path": 1}
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 2:
                rows.append((Path(f[0]).stem, f[0]))
                continue
            acc, p = f[cols["accession"]], f[cols["path"]]
            rows.append((acc, p))
            if "species" in cols and len(f) > cols["species"]:
                species[acc] = f[cols["species"]]
    return rows, species


ANCHOR_COLS = ["anchor_id", "genome", "species", "organism", "contig", "contig_len", "family", "model",
               "all_models", "all_families", "ambiguous", "score", "evalue",
               "anchor_strand", "anchor_len", "anchor_partial", "win_start", "win_end",
               "anchor_offset", "left_flank_bp", "right_flank_bp",
               "truncated_left", "truncated_right", "window_complete"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes", required=True, type=Path)
    ap.add_argument("--config", type=Path,
                    default=Path(__file__).resolve().parents[1] / "configs/blind/anchors.json")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--windows", choices=["none", "all"], default="none",
                    help="'none' (census, default) or 'all' (carve and store sequence)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--scratch", default=None)
    ap.add_argument("--limit", type=int, default=0)
    # A whole-corpus census is ~13 h on one node; sharding turns it into a job
    # array so a single failure costs one slice rather than the entire run.
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    args = ap.parse_args()
    require_slurm("anchor census")

    cfg = load_config(args.config)
    want_seq = args.windows == "all"
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "proteins").mkdir(exist_ok=True)

    anchor_hmm = args.out / "anchors.hmm"
    model_family = build_anchor_hmm(cfg, anchor_hmm)
    print(f"[census] {len(model_family)} models over {len(cfg['families'])} families",
          file=sys.stderr)

    genomes, species = read_genomes(args.genomes)
    if args.limit:
        genomes = genomes[:args.limit]
    if args.nshards > 1:
        total = len(genomes)
        # Stride rather than block: genome size varies systematically along the
        # manifest, so striding gives every shard a comparable workload.
        genomes = genomes[args.shard::args.nshards]
        print(f"[census] shard {args.shard}/{args.nshards}: "
              f"{len(genomes)}/{total} genomes", file=sys.stderr)
    print(f"[census] {len(genomes)} genomes, {args.workers} workers, windows={args.windows}",
          file=sys.stderr)

    fam_stats: dict[str, dict] = defaultdict(
        lambda: {"anchors": 0, "genomes": set(), "contigs": set(), "species": set(),
                 "complete_windows": 0, "partial_anchors": 0, "ambiguous": 0})
    failures: dict[str, str] = {}
    faa = {}
    n_anchor = 0

    at = open(args.out / "anchors.tsv", "w")
    at.write("\t".join(ANCHOR_COLS) + "\n")
    wj = open(args.out / "windows.jsonl", "w") if want_seq else None

    ctx = mp.get_context("fork")
    with ctx.Pool(args.workers, initializer=_init,
                  initargs=(cfg, str(anchor_hmm), model_family, args.scratch, want_seq)) as pool:
        for i, (acc, rows, windows, err) in enumerate(
                pool.imap_unordered(_one, genomes, chunksize=1), 1):
            if err:
                failures[acc] = err.strip().splitlines()[-1]
            sp = species.get(acc, "")
            # The manifest writes the literal "unknown" for genomes outside the
            # species views, which is truthy — treat it as absent.
            if sp in ("", "unknown", "NA"):
                sp = ""
            for r in rows:
                # View label when the genome is in one; otherwise the organism
                # read from the assembly's own header.
                r["species"] = sp or r.get("organism") or "unknown"
                at.write("\t".join(str(r[c]) for c in ANCHOR_COLS) + "\n")
                fam = r["family"]
                s = fam_stats[fam]
                s["anchors"] += 1
                s["genomes"].add(acc)
                s["contigs"].add(f"{acc}|{r['contig']}")
                # The resolved label, not the view label: `sp` is empty for the
                # 40% of genomes outside any view, which would collapse them all
                # to one pseudo-species and under-report diversity.
                s["species"].add(r["species"])
                s["complete_windows"] += int(r["window_complete"])
                s["partial_anchors"] += int(r["anchor_partial"])
                s["ambiguous"] += int(r["ambiguous"])
                if fam not in faa:
                    faa[fam] = open(args.out / "proteins" / f"{fam}.faa", "w")
                faa[fam].write(f">{r['anchor_id']} {fam} {r['model']} species={sp}\n"
                               f"{r['aa']}\n")
                n_anchor += 1
            if wj:
                for w in windows:
                    wj.write(json.dumps(w) + "\n")
            if i % 200 == 0 or i == len(genomes):
                print(f"[census] {i}/{len(genomes)} genomes, {n_anchor} anchors",
                      file=sys.stderr, flush=True)
    at.close()
    for fh in faa.values():
        fh.close()
    if wj:
        wj.close()

    with open(args.out / "summary_by_anchor_family.tsv", "w") as fh:
        fh.write("family\tanchors\tgenomes\tcontigs\tspecies\tcomplete_window_frac\t"
                 "partial_anchor_frac\tambiguous_frac\tanchors_per_genome\n")
        for fam, s in sorted(fam_stats.items(), key=lambda kv: -kv[1]["anchors"]):
            n = s["anchors"]
            fh.write(f"{fam}\t{n}\t{len(s['genomes'])}\t{len(s['contigs'])}\t"
                     f"{len(s['species'])}\t{s['complete_windows']/n:.4f}\t"
                     f"{s['partial_anchors']/n:.4f}\t{s['ambiguous']/n:.4f}\t"
                     f"{n/max(len(s['genomes']),1):.3f}\n")

    n_stub = sum(1 for v in failures.values() if v.startswith("empty-stub"))
    summary = {
        "n_genomes": len(genomes),
        "n_anchors": n_anchor,
        "windows_mode": args.windows,
        "per_family": {f: {"anchors": s["anchors"], "genomes": len(s["genomes"]),
                           "species": len(s["species"]),
                           "complete_window_frac": round(s["complete_windows"] / s["anchors"], 4)}
                       for f, s in fam_stats.items()},
        "n_empty_stubs": n_stub,
        "n_failures": len(failures) - n_stub,
        "failures": {k: v for k, v in failures.items() if not v.startswith("empty-stub")},
        "flank_bp": cfg["window"]["flank_bp"],
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary["per_family"], indent=2), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
