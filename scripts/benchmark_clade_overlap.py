#!/usr/bin/env python
"""Where do the known retrons land among the blind RT clades?

Benchmark branch — see data/known_retrons/README.md. This reads blind output
read-only and writes only under the benchmark directory. It must never be made
a dependency of `selected_clades.tsv` or of window generation.

    benchmark_clade_overlap.py --bench <benchmark_dir> --clades <clades_dir> \\
                               --out <benchmark_dir>

Answers three things the per-anchor table cannot:

  * how many of the blind-selected RT clades contain a known retron at all;
  * whether the known retrons concentrate in one clade or scatter across many —
    the interesting case is scatter, since a single clade would mean the blind
    decomposition had simply rediscovered "the retron cluster";
  * the known vs unknown split *inside* each selected clade, which is what a
    step-4 alignment would actually be built from.

A clade containing a known retron is not a success and a clade without one is
not a failure: this measures overlap, not recall of the discovery method. The
blind pipeline never used retron information to build these clades.
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        cols = next(fh).rstrip("\n").split("\t")
        return [dict(zip(cols, line.rstrip("\n").split("\t"))) for line in fh]


# --------------------------------------------------------------------------- #
# step-4 comparison
#
# Does a block the blind caller found sit where a known retron ncRNA actually
# is? The blocks were called with no access to RetronDB, so this is a scoring
# of frozen predictions, not a fit.

def load_window_frames(windows: Path) -> dict[str, tuple[str, int]]:
    """anchor_id -> (clade, anchor_offset), to put BLAST hits in the anchor frame.

    `annotate_known_retrons.py` reports hit coordinates in *window* space; blocks
    are in *anchor-frame* space (x = window_pos - anchor_offset). Without this
    shift the two are silently off by a few kb and nothing ever overlaps.
    """
    frames: dict[str, tuple[str, int]] = {}
    for p in sorted(windows.glob("*.fasta")):
        with open(p) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                f = line[1:].split()
                off = 0
                for kv in f[1:]:
                    if kv.startswith("anchor_offset="):
                        off = int(kv.split("=", 1)[1])
                frames[f[0]] = (p.stem, off)
    return frames


def overlap(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def compare_blocks(bench: list[dict], blocks_dir: Path, windows: Path,
                   clades_with_known: set[str], out: Path, slack: int) -> dict:
    blocks = read_tsv(blocks_dir / "blocks.tsv")
    members = read_tsv(blocks_dir / "block_members.tsv")
    frames = load_window_frames(windows)

    # Known ncRNA hits, shifted into the anchor frame of their own window.
    known_iv: dict[str, list[tuple[int, int, dict]]] = defaultdict(list)
    for r in bench:
        aid = r["anchor_id"]
        if aid not in frames or not r["hit_start_in_window"]:
            continue
        off = frames[aid][1]
        s, e = int(r["hit_start_in_window"]), int(r["hit_end_in_window"])
        known_iv[aid].append((min(s, e) - 1 - off, max(s, e) - off, r))

    by_block: dict[str, list[dict]] = defaultdict(list)
    for m in members:
        by_block[m["block_id"]].append(m)

    # Rank blocks within each clade by the blind composite score. Rank 1 is the
    # block the pipeline would hand a bench scientist first.
    # Seed with every real clade, including those that produced no block at all.
    # Reporting only block-bearing clades would silently shrink the denominator
    # and turn "found nothing here" into "not counted".
    per_clade: dict[str, list[dict]] = {p.stem: [] for p in windows.glob("*.fasta")}
    for b in blocks:
        if b["decoy"] == "True":
            continue
        per_clade.setdefault(b["clade"], []).append(b)
    for cl in per_clade:
        per_clade[cl].sort(key=lambda b: -float(b["composite_score"]))

    rows, clade_rows = [], []
    for cl, bs in sorted(per_clade.items()):
        first_rank, first_block, best_ov = None, "", 0
        n_hit_blocks = 0
        for rank, b in enumerate(bs, 1):
            hit_members, ov_bp = 0, 0
            for m in by_block.get(b["block_id"], []):
                for ks, ke, _ in known_iv.get(m["anchor_id"], []):
                    o = overlap(int(m["start"]) - slack, int(m["end"]) + slack, ks, ke)
                    if o > 0:
                        hit_members += 1
                        ov_bp = max(ov_bp, o)
                        break
            if hit_members:
                n_hit_blocks += 1
                if first_rank is None:
                    first_rank, first_block, best_ov = rank, b["block_id"], ov_bp
            rows.append({
                "clade": cl, "block_id": b["block_id"], "rank_in_clade": rank,
                "composite_score": b["composite_score"],
                "recurrence": b["recurrence"],
                "coord_tightness": b["coord_tightness"],
                "strand_consistency": b["strand_consistency"],
                "length_consistency": b["length_consistency"],
                "median_start": b["median_start"], "median_len": b["median_len"],
                "n_members": b["n_members"],
                "members_overlapping_known": hit_members,
                "max_overlap_bp": ov_bp,
                "overlaps_known_retron": bool(hit_members),
            })
        clade_rows.append({
            "clade": cl,
            "retron_positive": cl in clades_with_known,
            "n_blocks": len(bs),
            "n_blocks_overlapping_known": n_hit_blocks,
            "rank_of_first_overlapping_block": first_rank if first_rank else "",
            "first_overlapping_block": first_block,
            "max_overlap_bp": best_ov,
            "top_block_score": bs[0]["composite_score"] if bs else "",
        })

    bcols = ["clade", "block_id", "rank_in_clade", "composite_score", "recurrence",
             "coord_tightness", "strand_consistency", "length_consistency",
             "median_start", "median_len", "n_members",
             "members_overlapping_known", "max_overlap_bp", "overlaps_known_retron"]
    with open(out / "benchmark_block_overlap.tsv", "w") as fh:
        fh.write("\t".join(bcols) + "\n")
        for r in sorted(rows, key=lambda r: (r["clade"], r["rank_in_clade"])):
            fh.write("\t".join(str(r[c]) for c in bcols) + "\n")

    ccols = ["clade", "retron_positive", "n_blocks", "n_blocks_overlapping_known",
             "rank_of_first_overlapping_block", "first_overlapping_block",
             "max_overlap_bp", "top_block_score"]
    with open(out / "benchmark_block_rank.tsv", "w") as fh:
        fh.write("\t".join(ccols) + "\n")
        for r in sorted(clade_rows, key=lambda r: (not r["retron_positive"], r["clade"])):
            fh.write("\t".join(str(r[c]) for c in ccols) + "\n")

    def dist(bs: list[dict]) -> dict:
        sc = sorted(float(b["composite_score"]) for b in bs)
        if not sc:
            return {"n": 0}
        return {"n": len(sc), "median": round(sc[len(sc) // 2], 4),
                "p90": round(sc[int(len(sc) * 0.9)], 4), "max": round(sc[-1], 4)}

    rt_pos = [b for b in blocks if b["decoy"] == "False"
              and b["clade"] in clades_with_known]
    rt_unk = [b for b in blocks if b["decoy"] == "False" and b["family"] == "RT"
              and b["clade"] not in clades_with_known]
    decoy = [b for b in blocks if b["decoy"] == "True"]
    decoy_rt = [b for b in decoy if b["family"] == "RT"]

    pos = [r for r in clade_rows if r["retron_positive"]]
    ranks = [r["rank_of_first_overlapping_block"] for r in pos
             if r["rank_of_first_overlapping_block"]]
    return {
        "slack_bp": slack,
        "retron_positive_clades": len(pos),
        "retron_positive_clades_with_any_block": sum(1 for r in pos if r["n_blocks"]),
        "retron_positive_clades_with_no_block_at_all": sum(
            1 for r in pos if not r["n_blocks"]),
        "retron_positive_clades_with_an_overlapping_block": len(ranks),
        "rank_of_first_overlapping_block": sorted(ranks),
        "median_rank": statistics.median(ranks) if ranks else None,
        "rank_1_hits": sum(1 for r in ranks if r == 1),
        "score_distribution": {
            "retron_positive_RT_blocks": dist(rt_pos),
            "unknown_RT_blocks": dist(rt_unk),
            "decoy_RT_blocks": dist(decoy_rt),
            "all_decoy_blocks": dist(decoy),
        },
        "caveats": [
            "Blocks were called blind; this only scores frozen predictions.",
            "Overlap is to a blastn hit of a known ncRNA, not to a curated "
            "msr/msd boundary, so the comparison inherits RetronDB's limits.",
            "Rank is within a clade by blind composite score; a rank far from 1 "
            "means the pipeline found the retron but did not prioritise it.",
        ],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", required=True, type=Path)
    ap.add_argument("--clades", required=True, type=Path)
    ap.add_argument("--family", default="RT")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--blocks", type=Path,
                    help="step-4 output dir; enables the block-level comparison")
    ap.add_argument("--windows", type=Path,
                    help="window dir used for step 4, for the anchor-frame shift")
    ap.add_argument("--slack", type=int, default=100,
                    help="bp of tolerance when calling a block/ncRNA overlap")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    bench = read_tsv(args.bench / "benchmark_known_retrons.tsv")
    hit_anchor = {r["anchor_id"]: r for r in bench}
    print(f"{len(hit_anchor)} RT anchors matched a known retron")
    if not hit_anchor:
        print("nothing matched; clade overlap is empty by construction")

    # Selected clades of this family, at the recommended rung only — the same
    # file step 4 consumes, so the overlap describes the real step-4 input.
    sel = [r for r in read_tsv(args.clades / "selected_clades.tsv")
           if r["family"] == args.family]
    want = {(r["threshold"], r["cluster_id"]): r for r in sel}
    print(f"{len(sel)} selected {args.family} clades")

    members: dict[tuple, list[str]] = defaultdict(list)
    with open(args.clades / "cluster_members.tsv") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == args.family and (f[1], f[2]) in want:
                members[(f[1], f[2])].append(f[3])

    rows = []
    for key, mem in members.items():
        c = want[key]
        hits = [hit_anchor[a] for a in mem if a in hit_anchor]
        retrons = Counter(h["retron_id"] for h in hits if h["retron_id"])
        names = Counter(h["retron_name"] for h in hits if h["retron_name"])
        subtypes = Counter(h["retron_subtype"] for h in hits if h["retron_subtype"])
        rows.append({
            "family": args.family,
            "threshold": key[0],
            "cluster_id": key[1],
            "n_anchors": len(mem),
            "n_known_retron_anchors": len(hits),
            "n_unknown_anchors": len(mem) - len(hits),
            "known_frac": round(len(hits) / len(mem), 4) if mem else 0.0,
            "n_distinct_known_retrons": len(retrons),
            "retron_subtypes": ",".join(sorted(subtypes)) or "",
            "top_retron_names": ",".join(n for n, _ in names.most_common(3)),
            "median_identity": c.get("median_identity", ""),
            "n_complete_windows": c.get("n_complete_windows", ""),
            "rt_subtype": c.get("rt_subtype", ""),
        })

    cols = ["family", "threshold", "cluster_id", "n_anchors",
            "n_known_retron_anchors", "n_unknown_anchors", "known_frac",
            "n_distinct_known_retrons", "retron_subtypes", "top_retron_names",
            "median_identity", "n_complete_windows", "rt_subtype"]
    out_tsv = args.out / "benchmark_clade_overlap.tsv"
    with open(out_tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: -r["n_known_retron_anchors"]):
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    with_known = [r for r in rows if r["n_known_retron_anchors"]]
    in_selected = sum(r["n_known_retron_anchors"] for r in rows)
    all_retrons = {h["retron_id"] for h in hit_anchor.values() if h["retron_id"]}
    in_sel_retrons = set()
    for key, mem in members.items():
        for a in mem:
            if a in hit_anchor and hit_anchor[a]["retron_id"]:
                in_sel_retrons.add(hit_anchor[a]["retron_id"])

    summary = {
        "family": args.family,
        "matched_anchors_total": len(hit_anchor),
        "distinct_known_retrons_matched": len(all_retrons),
        "selected_clades": len(sel),
        "selected_clades_with_a_known_retron": len(with_known),
        "known_retron_anchors_inside_selected_clades": in_selected,
        "known_retron_anchors_outside_selected_clades": len(hit_anchor) - in_selected,
        "distinct_known_retrons_inside_selected_clades": len(in_sel_retrons),
        # One clade would mean the blind decomposition simply recovered "the
        # retron cluster"; several means retrons are spread across RT diversity.
        "concentration": ("none" if not with_known
                          else "single_clade" if len(with_known) == 1
                          else f"{len(with_known)}_clades"),
        "caveats": [
            "Overlap, not recall: a blastn hit shows a known ncRNA lies inside "
            "an RT window, not that blind discovery would have called it.",
            "Only clades at the recommended rung are counted, matching the "
            "selected_clades.tsv that step 4 consumes.",
            "Known-retron status was never used to build these clades.",
        ],
    }
    (args.out / "benchmark_clade_overlap.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: v for k, v in summary.items() if k != "caveats"}, indent=2))

    if with_known:
        w = ["cluster_id", "n_anchors", "n_known_retron_anchors",
             "n_unknown_anchors", "known_frac", "n_distinct_known_retrons",
             "median_identity", "retron_subtypes"]
        print()
        print("  ".join(c for c in w))
        for r in sorted(with_known, key=lambda r: -r["n_known_retron_anchors"]):
            print("  ".join(str(r[c])[:38] for c in w))
    if args.blocks and args.windows:
        print("\n=== step-4 block comparison (frozen blind predictions) ===")
        known_clades = {r["cluster_id"] for r in with_known}
        # Clade directory names encode family__thrX__cluster_id, so match on the
        # sanitized cluster id that extract_windows.py used for the filenames.
        safe = {_safe(cid) for cid in known_clades}
        pos_clades = {p.stem for p in args.windows.glob("*.fasta")
                      if p.stem.split("__", 2)[-1] in safe}
        bres = compare_blocks(bench, args.blocks, args.windows, pos_clades,
                              args.out, args.slack)
        (args.out / "benchmark_block_summary.json").write_text(
            json.dumps(bres, indent=2))
        print(json.dumps({k: v for k, v in bres.items() if k != "caveats"}, indent=2))
    print(f"-> {out_tsv}")
    return 0


def _safe(s: str) -> str:
    """Mirror extract_windows.py so clade ids map onto window filenames."""
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in s)[:80]


if __name__ == "__main__":
    raise SystemExit(main())
