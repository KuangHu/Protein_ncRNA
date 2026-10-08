#!/usr/bin/env python
"""Two shortlists over the triage tables: anchor-abutting blocks, and the
discovery top 10.

    discovery_shortlist.py --triage <triage_dir> --out <dir>

Neither is a new analysis. Both are views of `triage_blocks.tsv`, written down
so the selection rule is a file rather than a remembered filter.

`anchor_abutting.tsv` -- score >= 0.90, non-Cas3, and the block's near edge
sitting on the anchor boundary. A block edge that coincides with a start or stop
codon to within a few bases, with an interquartile range of ~0 across hundreds of
genomes, is a cis-module geometry: the element's RNA and its ORF share a
boundary. It is the signature the GroupII_RT blocks show, and it is strong enough
to be worth looking for on its own, independent of covariation power.

`discovery_top10.tsv` -- the same table with GroupII_RT, Cas3 and known-retron
overlaps removed, ranked. This is the entry point for deep annotation or
wet-lab: what is left after the families that are either positive-control-like
(GroupII_RT), flagged (Cas3), or already described (RetronDB) are set aside.

A note on `anchor_gap_bp`, because the obvious column is the wrong one. The
`distance_to_anchor` column from step 4 is the offset of the block's *start*:
for a downstream block 0 means it abuts the stop codon, but for an upstream
block it is the far edge, and 0 there is impossible by construction. Filtering
on `distance_to_anchor == 0` therefore silently selects downstream blocks only
and drops exactly the upstream blocks that motivate the filter -- including the
0.994 GroupII_RT pair, which ends at x=0 with IQR 0 and starts 549 bp away.
So this computes the gap at whichever edge faces the anchor:

    upstream:    gap = -relative_end_median        (0 = abuts the start codon)
    downstream:  gap = relative_start_median - anchor_len   (0 = abuts the stop)
    both/split:  the smaller of the two

with `anchor_edge_iqr_bp` giving the spread of that same edge across members,
since the tightness is half the signal.

Blind: coordinates, ORF geometry and each block's own fold. Benchmark columns
are read only to exclude already-described loci from the discovery list, which
is a filter on the output and never on the score.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

# How close an edge has to sit to count as abutting. Not zero: the block
# boundary is an alignment edge and the ORF boundary is a Prodigal call, so a
# few bases of disagreement are expected even when they are the same boundary.
ABUT_BP = 10
# And how tight that edge has to be across members. A single genome can place an
# edge anywhere; the claim is about the clade.
ABUT_IQR_BP = 25

VERDICT_RANK = {"covariation_supported": 0, "covariation_borderline": 1,
                "tested_no_covariation": 2, "underpowered": 3,
                "covariation_test_failed": 4}


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        cols = next(fh).rstrip("\n").split("\t")
        return [dict(zip(cols, line.rstrip("\n").split("\t"))) for line in fh]


def write_tsv(path: Path, rows: list[dict], cols: list[str]) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def num(r: dict, k: str, default=0.0) -> float:
    try:
        return float(r[k])
    except (KeyError, ValueError, TypeError):
        return default


def anchor_gap(r: dict) -> tuple[int, int, str]:
    """(gap_bp, edge_iqr_bp, which_edge) at the edge facing the anchor."""
    alen = num(r, "anchor_len")
    up_gap = -num(r, "relative_end_median")
    up_iqr = abs(num(r, "relative_end_q75") - num(r, "relative_end_q25"))
    dn_gap = num(r, "relative_start_median") - alen
    dn_iqr = abs(num(r, "relative_start_q75") - num(r, "relative_start_q25"))
    side = r.get("side_of_anchor", "")
    if side == "upstream":
        return int(up_gap), int(up_iqr), "end_at_start_codon"
    if side == "downstream":
        return int(dn_gap), int(dn_iqr), "start_at_stop_codon"
    # flanks_both_sides / split: report whichever edge is nearer the anchor.
    if abs(up_gap) <= abs(dn_gap):
        return int(up_gap), int(up_iqr), "end_at_start_codon"
    return int(dn_gap), int(dn_iqr), "start_at_stop_codon"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--triage", required=True, type=Path,
                    help="triage_report.py output dir")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--min-score", type=float, default=0.90)
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--exclude-families", nargs="*",
                    default=["GroupII_RT", "Cas3"])
    args = ap.parse_args()
    out = args.out or args.triage
    out.mkdir(parents=True, exist_ok=True)

    rows = read_tsv(args.triage / "triage_blocks.tsv")
    for r in rows:
        g, iqr, edge = anchor_gap(r)
        r["anchor_gap_bp"], r["anchor_edge_iqr_bp"], r["abutting_edge"] = g, iqr, edge
        r["abuts_anchor"] = abs(g) <= ABUT_BP and iqr <= ABUT_IQR_BP
    print(f"{len(rows)} blocks read")

    # ---- shortlist 1: high-score anchor-abutting -------------------------
    ab = [r for r in rows
          if num(r, "score") >= args.min_score
          and r["family"] != "Cas3"
          and r["abuts_anchor"]]
    ab.sort(key=lambda r: -num(r, "score"))
    cols1 = ["family", "clade_id", "block_id", "score", "side_of_anchor",
             "abutting_edge", "anchor_gap_bp", "anchor_edge_iqr_bp",
             "relative_start_median", "relative_end_median", "anchor_len",
             "n_members", "member_fraction", "covariation_verdict",
             "rscape_power", "TP_covarying_pairs", "FP_pairs",
             "padded_coding_frac", "median_pairwise_identity",
             "same_clade_neighbor_blocks_500bp", "nearby_block_group_id",
             "benchmark_overlaps_known_ncrna"]
    write_tsv(out / "anchor_abutting.tsv", ab, cols1)
    print(f"anchor-abutting (score>={args.min_score}, non-Cas3, |gap|<={ABUT_BP}, "
          f"IQR<={ABUT_IQR_BP}): {len(ab)} blocks")

    # ---- shortlist 2: discovery top N ------------------------------------
    def known(r):
        return str(r.get("benchmark_overlaps_known_ncrna", "")).lower() in (
            "true", "1", "yes")

    pool = [r for r in rows
            if r["family"] not in args.exclude_families and not known(r)]
    # Rank on the blind columns, in the order they carry weight: covariation
    # first where it exists, then score, then how much of the clade carries the
    # block, then how close it sits to the anchor, then how noncoding it is.
    #
    # `covariation_borderline` is not one thing, so it cannot be one rank. The
    # verdict only asks whether TP > 0; a block with TP=1 and FP=118 (PPV 0.84%)
    # earns the same label as one with TP=4 and FP=0, and sorting on score alone
    # puts the first above the second. Low-PPV borderline blocks are noise with a
    # covarying pair in it, so they sort after clean ones of the same verdict
    # rather than being interleaved with them by score.
    for r in rows:
        r["_ppv_band"] = (0 if r["covariation_verdict"] not in
                          ("covariation_borderline", "covariation_supported")
                          else 0 if num(r, "PPV", 0.0) >= 50.0 else 1)
    pool.sort(key=lambda r: (
        VERDICT_RANK.get(r["covariation_verdict"], 9),
        r["_ppv_band"],
        -num(r, "score"),
        -num(r, "member_fraction"),
        abs(r["anchor_gap_bp"]),
        num(r, "padded_coding_frac"),
    ))
    top = pool[:args.top]
    cols2 = ["rank", "family", "clade_id", "block_id", "score",
             "covariation_verdict", "rscape_power", "rscape_expected_covarying",
             "TP_covarying_pairs", "FP_pairs", "PPV", "member_fraction",
             "n_members", "side_of_anchor", "anchor_gap_bp",
             "anchor_edge_iqr_bp", "relative_start_median",
             "relative_end_median", "anchor_len", "padded_coding_frac",
             "median_pairwise_identity", "nearby_block_group_id",
             "same_clade_neighbor_blocks_500bp"]
    for i, r in enumerate(top, 1):
        r["rank"] = i
    write_tsv(out / "discovery_top10.tsv", top, cols2)
    print(f"discovery pool (excluding {', '.join(args.exclude_families)} and "
          f"known-retron overlaps): {len(pool)} blocks, top {len(top)} written")

    summary = {
        "n_blocks": len(rows),
        "abut_bp": ABUT_BP, "abut_iqr_bp": ABUT_IQR_BP,
        "min_score_anchor_abutting": args.min_score,
        "n_anchor_abutting": len(ab),
        "anchor_abutting_by_family": {
            f: sum(1 for r in ab if r["family"] == f)
            for f in sorted({r["family"] for r in ab})},
        "n_blocks_abutting_any_score": sum(1 for r in rows if r["abuts_anchor"]),
        "excluded_families": args.exclude_families,
        "n_known_retron_overlaps_excluded": sum(1 for r in rows if known(r)),
        "discovery_pool": len(pool),
        "discovery_pool_by_family": {
            f: sum(1 for r in pool if r["family"] == f)
            for f in sorted({r["family"] for r in pool})},
    }
    (out / "discovery_shortlist.json").write_text(json.dumps(summary, indent=2))
    print("\n" + json.dumps(summary, indent=2))
    print(f"-> {out}/anchor_abutting.tsv, {out}/discovery_top10.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
