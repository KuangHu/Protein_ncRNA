#!/usr/bin/env python
"""Split RT anchors into a known-retron benchmark set and a blind-unknown set.

Benchmark branch — see data/known_retrons/README.md. This must never run as part
of producing `selected_clades.tsv`; it reads blind output and partitions it
*afterwards*.

    annotate_known_retrons.py --known data/known_retrons/known_retrons.tsv \\
                              --census <corrected> --windows out/rt_windows \\
                              --out out/known_retron_benchmark

Writes `benchmark_known_retrons.tsv` (RT anchors overlapping a known retron
ncRNA) and `blind_unknown_rt_anchors.tsv` (everything else — the honest blind
discovery set).

Two matching routes, with very different yields:

  * **accession** — attempted and expected to recover nothing. RetronDB is keyed
    by PATRIC or RefSeq protein accessions; this pipeline calls ORFs de novo and
    `genomes_db` is keyed by assembly accession, so there is no join. Reported
    anyway so the gap stays visible.
  * **sequence** — blastn of the known ncRNAs against extracted RT windows. This
    is the route that works, against **180** sequence-bearing rows out of 1,928.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_ncrna.compute import require_slurm
from protein_ncrna.tools import which


def load_known(path: Path) -> tuple[list[dict], dict[str, dict]]:
    rows, by_acc = [], {}
    with open(path) as fh:
        c = {k: i for i, k in enumerate(next(fh).rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            r = {k: f[i] for k, i in c.items()}
            rows.append(r)
            if r["accession"]:
                by_acc[r["accession"]] = r
    return rows, by_acc


def run_blast(query: Path, subject_fastas: list[Path], out: Path, threads: int,
              evalue: float, min_ident: float, min_cov: float,
              known_len: dict[str, int]) -> dict[str, list[dict]]:
    """blastn known ncRNAs against window FASTAs. Returns anchor_id -> hits."""
    blastn = which("blastn")
    makedb = which("makeblastdb") if _has("makeblastdb") else None
    hits: dict[str, list[dict]] = defaultdict(list)

    for sf in subject_fastas:
        if sf.stat().st_size == 0:
            continue
        res = out / f"blast_{sf.stem}.tsv"
        cmd = [blastn, "-query", str(query), "-subject", str(sf),
               "-outfmt", "6 qseqid sseqid pident length qlen sstart send evalue bitscore",
               "-evalue", str(evalue), "-num_threads", str(threads),
               # Known retron ncRNAs are short and structured; the default
               # megablast word size misses diverged homologues entirely.
               "-task", "blastn", "-word_size", "11"]
        if makedb:
            db = out / f"db_{sf.stem}"
            subprocess.run([makedb, "-in", str(sf), "-dbtype", "nucl", "-out", str(db)],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            cmd = [blastn, "-query", str(query), "-db", str(db),
                   "-outfmt", "6 qseqid sseqid pident length qlen sstart send evalue bitscore",
                   "-evalue", str(evalue), "-num_threads", str(threads),
                   "-task", "blastn", "-word_size", "11"]
        proc = subprocess.run(cmd, stdout=open(res, "w"), stderr=subprocess.PIPE)
        if proc.returncode != 0:
            print(f"  blastn failed on {sf.name}: "
                  f"{proc.stderr.decode(errors='replace')[-200:]}", file=sys.stderr)
            continue
        with open(res) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                q, s, pid, aln = f[0], f[1], float(f[2]), int(f[3])
                qlen = known_len.get(q, int(f[4]) if f[4].isdigit() else 0)
                cov = aln / qlen if qlen else 0.0
                if pid < min_ident or cov < min_cov:
                    continue
                hits[s].append({"retron_id": q, "pident": pid, "aln_len": aln,
                                "qcov": round(cov, 3),
                                "sstart": int(f[5]), "send": int(f[6]),
                                "evalue": float(f[7]), "bitscore": float(f[8])})
    return hits


def _has(name: str) -> bool:
    try:
        which(name)
        return True
    except FileNotFoundError:
        return False


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--known", required=True, type=Path)
    ap.add_argument("--census", required=True, type=Path)
    ap.add_argument("--windows", type=Path,
                    help="directory of window FASTAs from extract_windows.py")
    ap.add_argument("--family", default="RT")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--min-ident", type=float, default=70.0)
    ap.add_argument("--min-cov", type=float, default=0.50)
    args = ap.parse_args()
    require_slurm("known-retron benchmark")
    args.out.mkdir(parents=True, exist_ok=True)

    known, by_acc = load_known(args.known)
    with_seq = [r for r in known if r["ncrna"]]
    known_len = {r["retron_id"]: int(r["ncrna_len"]) for r in known}
    print(f"{len(known)} known retrons, {len(with_seq)} with a sequence")

    anchors = []
    with open(args.census / "anchors.tsv") as fh:
        hdr = next(fh).rstrip("\n").split("\t")
        c = {k: i for i, k in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[c["family"]] == args.family:
                anchors.append({k: f[i] for k, i in c.items()})
    print(f"{len(anchors)} {args.family} anchors")

    # Route 1: accession. Expected to recover nothing; measured, not assumed.
    acc_hits = {a["anchor_id"]: by_acc[a["genome"]] for a in anchors
                if a["genome"] in by_acc}
    print(f"accession matches: {len(acc_hits)} "
          f"(expected ~0 — RetronDB uses PATRIC/protein ids)")

    # Route 2: sequence.
    seq_hits: dict[str, list[dict]] = {}
    if args.windows:
        fastas = sorted(args.windows.glob("*.fasta"))
        query = args.known.with_name("known_retron_ncrna.fasta")
        if not query.exists():
            raise SystemExit(f"missing {query}; run import_retrondb.py first")
        print(f"blastn {len(with_seq)} ncRNAs against {len(fastas)} window FASTAs")
        seq_hits = run_blast(query, fastas, args.out, args.threads, args.evalue,
                             args.min_ident, args.min_cov, known_len)
        print(f"sequence matches: {len(seq_hits)} anchors")
    else:
        print("no --windows given; sequence matching skipped")

    matched = set(acc_hits) | set(seq_hits)
    by_id = {r["retron_id"]: r for r in known}

    bench_cols = ["anchor_id", "genome", "species", "family", "match_route",
                  "retron_id", "retron_name", "rt_clade", "retron_subtype",
                  "msr_msd_family", "pident", "qcov", "evalue",
                  "hit_start_in_window", "hit_end_in_window", "anchor_offset"]
    with open(args.out / "benchmark_known_retrons.tsv", "w") as fh:
        fh.write("\t".join(bench_cols) + "\n")
        for a in anchors:
            aid = a["anchor_id"]
            if aid not in matched:
                continue
            best = max(seq_hits.get(aid, []), key=lambda h: h["bitscore"], default=None)
            k = by_id.get(best["retron_id"]) if best else acc_hits.get(aid, {})
            fh.write("\t".join(str(x) for x in [
                aid, a["genome"], a["species"], a["family"],
                "sequence" if best else "accession",
                (k or {}).get("retron_id", ""), (k or {}).get("retron_name", ""),
                (k or {}).get("rt_clade", ""), (k or {}).get("retron_subtype", ""),
                (k or {}).get("msr_msd_family", ""),
                best["pident"] if best else "", best["qcov"] if best else "",
                best["evalue"] if best else "",
                best["sstart"] if best else "", best["send"] if best else "",
                a["anchor_offset"],
            ]) + "\n")

    with open(args.out / "blind_unknown_rt_anchors.tsv", "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for a in anchors:
            if a["anchor_id"] not in matched:
                fh.write("\t".join(a[k] for k in hdr) + "\n")

    recovered = {h["retron_id"] for hs in seq_hits.values() for h in hs}
    report = {
        "family": args.family,
        "anchors": len(anchors),
        "known_retrons": len(known),
        "known_with_sequence": len(with_seq),
        "matched_anchors": len(matched),
        "matched_by_accession": len(acc_hits),
        "matched_by_sequence": len(seq_hits),
        "distinct_known_retrons_recovered": len(recovered),
        "recall_vs_sequence_bearing": (round(len(recovered) / len(with_seq), 4)
                                       if with_seq else None),
        "blind_unknown_anchors": len(anchors) - len(matched),
        "caveats": [
            "Recall is against the %d sequence-bearing rows, not %d retrons."
            % (len(with_seq), len(known)),
            "A blastn hit shows a known ncRNA lies inside an RT window; it does "
            "not show the blind pipeline would have called that block.",
            "RetronDB and genomes_db are both enterobacteria-heavy, so agreement "
            "between them is not evidence of generality.",
        ],
        "thresholds": {"evalue": args.evalue, "min_ident": args.min_ident,
                       "min_cov": args.min_cov},
    }
    (args.out / "benchmark_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps({k: v for k, v in report.items() if k != "caveats"}, indent=2))
    print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
