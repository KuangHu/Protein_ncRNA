#!/usr/bin/env python
"""Final step-4 candidate table: blocks that clear the frozen score floor.

    candidate_table.py --blocks <k3_dir> --out <dir> [--benchmark <bench_dir>]

The table itself is blind. `--benchmark` adds columns prefixed `benchmark_`,
which are annotations on an already-frozen ranking: they are attached after the
score and the ordering are fixed, and no benchmark value is an input to any
other column. Without `--benchmark` the table is identical minus those columns.

`decoy_empirical_p_family` is the empirical FDR from the within-family
permutation null, per family, at this block's own score:

    decoy blocks in family with score >= s, divided by the decoy replicate count
    ----------------------------------------------------------------------------
    real blocks in family with score >= s

so 0.0 means no decoy clade in that family ever produced a block scoring this
well. It is a rate ratio against a matched null, not a parametric p-value, and
it is bounded below by the number of decoy clades: with 129 decoys, 0.0 means
"below the resolution of this null", not "impossible".
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

COLS = ["block_id", "family", "clade_id", "score", "high_confidence",
        "benchmark_only", "n_members", "member_fraction", "block_len",
        "relative_start", "relative_end", "n_references_contributing",
        "frac_members_with_noncoding_hsp", "decoy_empirical_p_family",
        "nearest_anchor_distance", "fasta_path"]
BENCH_COLS = ["benchmark_overlaps_known_ncrna", "benchmark_rank_in_clade",
              "benchmark_members_overlapping", "benchmark_max_overlap_bp"]


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        cols = next(fh).rstrip("\n").split("\t")
        return [dict(zip(cols, line.rstrip("\n").split("\t"))) for line in fh]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--blocks", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--config", type=Path,
                    default=Path(__file__).resolve().parents[1]
                    / "configs" / "blind" / "step4.json")
    ap.add_argument("--benchmark", type=Path,
                    help="benchmark dir; adds benchmark_* annotation columns only")
    ap.add_argument("--decoys", type=Path,
                    help="decoy window dir; reads decoy_manifest.json to label "
                         "each family's null strength")
    ap.add_argument("--allow-incomplete", action="store_true",
                    help="build the table from a step-4 run that lost clades; "
                         "the decoy denominator will be short and the empirical "
                         "FDR correspondingly optimistic")
    ap.add_argument("--replicates", type=int, default=3,
                    help="decoy replicates per real clade, for null scaling")
    args = ap.parse_args()
    cfg = json.loads(args.config.read_text())
    floor = cfg["score_floor"]
    hi = cfg["high_confidence_floor"]
    bench_only = set(cfg.get("benchmark_only_families", []))
    args.out.mkdir(parents=True, exist_ok=True)

    # Every number in this table is read against the decoy null, so a step-4 run
    # that lost clades cannot be used: a missing decoy clade shrinks the null's
    # denominator and lowers the empirical FDR for free. Refuse rather than
    # annotate, because the resulting table would look entirely normal.
    ds = args.blocks / "discovery_summary.json"
    if ds.exists():
        d = json.loads(ds.read_text())
        if d.get("complete") is False and args.allow_incomplete:
            print(f"WARNING: {args.blocks} lost {d.get('n_failures')} clades; "
                  f"the decoy denominator is short and every "
                  f"decoy_empirical_p_family below is optimistic",
                  file=sys.stderr)
        elif d.get("complete") is False:
            raise SystemExit(
                f"{args.blocks} is an incomplete step-4 run "
                f"({d.get('n_failures')} clades failed; see failure_report.tsv). "
                f"Re-run step 4, or pass --allow-incomplete to build a table "
                f"whose decoy denominator is known to be short.")
        if "complete" not in d:
            print(f"note: {ds} predates the completeness flag; "
                  f"clade counts not verified", file=sys.stderr)
    elif not args.allow_incomplete:
        raise SystemExit(f"no {ds}; cannot verify the step-4 run was complete")

    blocks = read_tsv(args.blocks / "blocks.tsv")
    stats = {s["clade"]: s for s in read_tsv(args.blocks / "clade_stats.tsv")}
    real = [b for b in blocks if b["decoy"] == "False"]
    decoy = [b for b in blocks if b["decoy"] == "True"]
    S = lambda b: float(b["composite_score"])

    def emp_fdr(fam: str, s: float) -> float:
        d = sum(1 for b in decoy if b["family"] == fam and S(b) >= s) / args.replicates
        r = sum(1 for b in real if b["family"] == fam and S(b) >= s)
        return round(d / r, 4) if r else ""

    # Benchmark annotations, keyed by block, attached after scoring.
    bann: dict[str, dict] = {}
    if args.benchmark:
        bp = args.benchmark / "benchmark_block_overlap.tsv"
        if bp.exists():
            for r in read_tsv(bp):
                bann[r["block_id"]] = r
        else:
            print(f"no {bp}; benchmark columns will be blank", file=sys.stderr)

    rows = []
    for b in real:
        s = S(b)
        if s < floor:
            continue
        cl = b["clade"]
        st = stats.get(cl, {})
        start = int(b["median_start"])
        blen = int(b["median_len"])
        r = {
            "block_id": b["block_id"],
            "family": b["family"],
            "clade_id": cl,
            "score": round(s, 4),
            "high_confidence": s >= hi,
            "benchmark_only": b["family"] in bench_only,
            "n_members": b["n_members"],
            "member_fraction": b.get("member_fraction", b["recurrence"]),
            "block_len": blen,
            # Anchor frame: 0 is the anchor's first base, negative is upstream.
            "relative_start": start,
            "relative_end": start + blen,
            "n_references_contributing": b.get("n_references_contributing", ""),
            "frac_members_with_noncoding_hsp":
                st.get("frac_members_with_noncoding_hsp", ""),
            "decoy_empirical_p_family": emp_fdr(b["family"], s),
            # Signed distance to the anchor ORF: negative upstream of the start
            # codon, positive past the stop, 0 if the block abuts or overlaps.
            "nearest_anchor_distance": b["dist_to_anchor"],
            "fasta_path": str(args.blocks / "blocks" / f"{b['block_id']}.fasta"),
        }
        if args.benchmark:
            a = bann.get(b["block_id"], {})
            r["benchmark_overlaps_known_ncrna"] = a.get("overlaps_known_retron", "")
            r["benchmark_rank_in_clade"] = a.get("rank_in_clade", "")
            r["benchmark_members_overlapping"] = a.get("members_overlapping_known", "")
            r["benchmark_max_overlap_bp"] = a.get("max_overlap_bp", "")
        rows.append(r)

    cols = COLS + (BENCH_COLS if args.benchmark else [])
    rows.sort(key=lambda r: (r["benchmark_only"], -r["score"]))
    out = args.out / "step4_candidates.tsv"
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    # Per-family decoy counts at both floors -- the null these rows are read
    # against, carried alongside them so the table is never quoted bare.
    # A family whose decoys could not be drawn from enough distinct source
    # clades has a null that retains relatedness. That biases its FDR *upward*,
    # so the column is an upper bound rather than a measurement, and it must be
    # labelled instead of quoted bare.
    null_strength: dict[str, str] = {}
    dm = args.decoys / "decoy_manifest.json" if args.decoys else None
    if dm and dm.exists():
        d = json.loads(dm.read_text())
        null_strength = {f: v.get("null_strength", "")
                         for f, v in d.get("per_family", {}).items()}
        if not null_strength and d.get("decoys"):
            print(f"note: {dm} predates null-strength diagnostics", file=sys.stderr)

    fams = sorted({b["family"] for b in real})
    fcols = ["family", "real_blocks_ge_080", "decoy_blocks_ge_080",
             "decoy_scaled_ge_080", "empirical_fdr_ge_080",
             "real_blocks_ge_090", "decoy_blocks_ge_090",
             "decoy_scaled_ge_090", "empirical_fdr_ge_090", "null_strength",
             "benchmark_only"]
    frows = []
    for fam in fams:
        row = {"family": fam, "benchmark_only": fam in bench_only,
               "null_strength": null_strength.get(fam, "unknown")}
        for tag, c in (("080", floor), ("090", hi)):
            r_n = sum(1 for b in real if b["family"] == fam and S(b) >= c)
            d_n = sum(1 for b in decoy if b["family"] == fam and S(b) >= c)
            row[f"real_blocks_ge_{tag}"] = r_n
            row[f"decoy_blocks_ge_{tag}"] = d_n
            row[f"decoy_scaled_ge_{tag}"] = round(d_n / args.replicates, 2)
            row[f"empirical_fdr_ge_{tag}"] = (round(d_n / args.replicates / r_n, 4)
                                              if r_n else "")
        frows.append(row)
    with open(args.out / "step4_family_null.tsv", "w") as fh:
        fh.write("\t".join(fcols) + "\n")
        for r in frows:
            fh.write("\t".join(str(r[c]) for c in fcols) + "\n")

    summary = {
        "frozen_config": {k: cfg[k] for k in
                          ("reference_picker", "k_references", "min_recurrence",
                           "score_floor", "high_confidence_floor")},
        "candidates_ge_score_floor": len(rows),
        "high_confidence": sum(1 for r in rows if r["high_confidence"]),
        "benchmark_only_rows": sum(1 for r in rows if r["benchmark_only"]),
        "discovery_candidates": sum(1 for r in rows if not r["benchmark_only"]),
        "per_family": frows,
        "benchmark_annotated": bool(args.benchmark),
        "families_with_degenerate_null": sorted(
            f for f, v in null_strength.items() if v == "degenerate"),
    }
    (args.out / "step4_candidates.json").write_text(json.dumps(summary, indent=2))

    w = ["family", "real_blocks_ge_080", "decoy_scaled_ge_080",
         "empirical_fdr_ge_080", "real_blocks_ge_090", "decoy_scaled_ge_090",
         "null_strength", "benchmark_only"]
    print("  ".join(c for c in w))
    for r in frows:
        print("  ".join(str(r[c]).ljust(len(c)) for c in w))
    print(f"\n{len(rows)} candidates >= {floor} "
          f"({summary['high_confidence']} >= {hi}, "
          f"{summary['benchmark_only_rows']} flagged benchmark_only)")
    print(f"-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
