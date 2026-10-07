#!/usr/bin/env python
"""Merge the per-shard outputs of an array census into one census directory.

    merge_shards.py --census /path/to/census_<arrayjobid>

Reads `<census>/shard_*/` and writes `anchors.tsv`, `proteins/<family>.faa`,
`summary_by_anchor_family.tsv` and `summary.json` at the top level, in the same
layout a single-shard run produces, so census_report.py cannot tell the
difference.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", required=True, type=Path)
    args = ap.parse_args()

    shards = sorted(args.census.glob("shard_*"))
    if not shards:
        print(f"no shard_* directories under {args.census}", file=sys.stderr)
        return 1
    print(f"merging {len(shards)} shards")

    (args.census / "proteins").mkdir(exist_ok=True)
    header = None
    n_rows = 0
    fam_stats: dict[str, dict] = defaultdict(
        lambda: {"anchors": 0, "genomes": set(), "contigs": set(), "species": set(),
                 "complete_windows": 0, "partial_anchors": 0, "ambiguous": 0})
    totals = {"n_genomes": 0, "n_empty_stubs": 0, "n_failures": 0}
    failures: dict[str, str] = {}
    missing = []

    with open(args.census / "anchors.tsv", "w") as out:
        for sh in shards:
            at = sh / "anchors.tsv"
            if not at.exists():
                missing.append(sh.name)
                continue
            with open(at) as fh:
                hdr = fh.readline().rstrip("\n")
                if header is None:
                    header = hdr
                    out.write(hdr + "\n")
                elif hdr != header:
                    print(f"{sh.name}: header mismatch — shards from different runs?",
                          file=sys.stderr)
                    return 1
                cols = {c: i for i, c in enumerate(header.split("\t"))}
                for line in fh:
                    out.write(line)
                    f = line.rstrip("\n").split("\t")
                    fam = f[cols["family"]]
                    s = fam_stats[fam]
                    s["anchors"] += 1
                    s["genomes"].add(f[cols["genome"]])
                    s["contigs"].add(f"{f[cols['genome']]}|{f[cols['contig']]}")
                    s["species"].add(f[cols["species"]])
                    s["complete_windows"] += f[cols["window_complete"]] == "True"
                    s["partial_anchors"] += f[cols["anchor_partial"]] == "True"
                    s["ambiguous"] += f[cols["ambiguous"]] == "True"
                    n_rows += 1

            sj = sh / "summary.json"
            if sj.exists():
                d = json.loads(sj.read_text())
                totals["n_genomes"] += d.get("n_genomes", 0)
                totals["n_empty_stubs"] += d.get("n_empty_stubs", 0)
                totals["n_failures"] += d.get("n_failures", 0)
                failures.update(d.get("failures", {}))

    if missing:
        print(f"WARNING: {len(missing)} shards produced no anchors.tsv: "
              f"{', '.join(missing)}", file=sys.stderr)

    # Per-family FASTA, concatenated across shards.
    for fam in fam_stats:
        with open(args.census / "proteins" / f"{fam}.faa", "w") as out:
            for sh in shards:
                src = sh / "proteins" / f"{fam}.faa"
                if src.exists():
                    out.write(src.read_text())

    with open(args.census / "summary_by_anchor_family.tsv", "w") as fh:
        fh.write("family\tanchors\tgenomes\tcontigs\tspecies\tcomplete_window_frac\t"
                 "partial_anchor_frac\tambiguous_frac\tanchors_per_genome\n")
        for fam, s in sorted(fam_stats.items(), key=lambda kv: -kv[1]["anchors"]):
            n = s["anchors"]
            fh.write(f"{fam}\t{n}\t{len(s['genomes'])}\t{len(s['contigs'])}\t"
                     f"{len(s['species'])}\t{s['complete_windows']/n:.4f}\t"
                     f"{s['partial_anchors']/n:.4f}\t{s['ambiguous']/n:.4f}\t"
                     f"{n/max(len(s['genomes']),1):.3f}\n")

    summary = {
        "n_shards": len(shards),
        "shards_missing_output": missing,
        "n_anchors": n_rows,
        **totals,
        "per_family": {f: {"anchors": s["anchors"], "genomes": len(s["genomes"]),
                           "species": len(s["species"]),
                           "complete_window_frac": round(s["complete_windows"]/s["anchors"], 4)}
                       for f, s in fam_stats.items()},
        "failures": dict(list(failures.items())[:50]),
    }
    (args.census / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"{n_rows} anchors over {totals['n_genomes']} genomes, "
          f"{len(fam_stats)} families -> {args.census}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
