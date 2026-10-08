#!/usr/bin/env python
"""Report the stratified Arm 2 run as four tables, one per lane's question.

    arm2_lane_report.py --blocks <blocks dir> --decoys <decoys dir> \\
                        --lanes <run_bag_lanes.tsv> --strata <bag_strata.tsv> \\
                        --floor <score_floor.json> --out <dir>

The lanes answer different questions and must not share a ranking:

  positive_control_recovery.tsv  IS110 / group II / IS200-IS605 / TnpB-Cas12.
      These are spike-ins. The number is RECOVERY: what fraction of insertion
      systems already known to carry a structured RNA does the blind caller
      find? They are calibration and are excluded from novelty ranking, so a
      high score here is reassurance, never a result.
  discovery_candidates.tsv       everything else with a block. The result.
  host_like_background.tsv       chemotaxis receptors, ushers, surface proteins.
      The false-positive pressure estimate. If these score like the
      transposases, the caller is broken and nothing else in here means
      anything.
  all_bags_annotated.tsv         every bag run, whatever it did.

`decoy_clean` is per bag and per null: the bag's own decoys stayed below the
floor. It is not a p-value and not an FDR -- with a handful of decoy clades per
bag it is an empirical upper bound, "no observed decoy support", and the decoy
count backing it is carried in the same row so it cannot be read as more.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

STRICT = ("within_bag_shuffle", "same_family_permutation")
MODES = ("genomic_background",) + STRICT

# Pfam novelty for the final sort key: higher is more interesting. A bag with no
# specific domain is the most novel thing the run can return; a bag that is a
# well-characterised non-RNA mobile element is the least.
NOVELTY = {"unknown": 3, "discovery": 2, "host_like": 1, "known_rna_mge": 0}


def mode_of(clade: str) -> str | None:
    for m in MODES:
        if clade.startswith(f"DECOY_{m}__"):
            return m
    return None


def bag_of(clade: str, bags: set[str]) -> str | None:
    for b in bags:
        if clade == b or clade.endswith(f"__{b}"):
            return b
    return None


def write(path: Path, cols: list[str], rows: list[dict]) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--blocks", required=True, type=Path)
    ap.add_argument("--decoys", required=True, type=Path)
    ap.add_argument("--lanes", required=True, type=Path)
    ap.add_argument("--strata", required=True, type=Path)
    ap.add_argument("--floor", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    fl = json.loads(args.floor.read_text())
    floor = fl.get("floor_no_decoy_support")
    floor_n = fl.get("floor_backed_by_n_decoy_blocks", 0)
    if floor is None:
        print("WARNING no floor could be calibrated; decoy_clean columns will "
              "be blank and nothing is reported as clean")

    lanes = {r["bag"]: r for r in csv.DictReader(
        open(args.lanes), delimiter="\t")}
    strata = {r["bag"]: r for r in csv.DictReader(
        open(args.strata), delimiter="\t")}
    bags = set(lanes)

    # Which bags step 4 actually scored, so a bag that failed extraction is
    # never reported as a bag that produced no block.
    scored: set[str] = set()
    cs = args.blocks / "clade_stats.tsv"
    if cs.exists():
        for r in csv.DictReader(open(cs), delimiter="\t"):
            if r.get("decoy") == "False":
                b = bag_of(r["clade"], bags)
                if b:
                    scored.add(b)

    built: dict[tuple[str, str], set[str]] = defaultdict(set)
    for p in sorted(args.decoys.glob("DECOY_*.fasta")):
        m, b = mode_of(p.stem), bag_of(p.stem, bags)
        if m and b:
            built[(b, m)].add(p.stem)

    real: dict[str, list[dict]] = defaultdict(list)
    dec: dict[tuple[str, str], list[float]] = defaultdict(list)
    dec_clades: dict[tuple[str, str], set[str]] = defaultdict(set)
    for r in csv.DictReader(open(args.blocks / "blocks.tsv"), delimiter="\t"):
        s = float(r["composite_score"])
        if r["decoy"] == "True":
            m, b = mode_of(r["clade"]), bag_of(r["clade"], bags)
            if m and b:
                dec[(b, m)].append(s)
                dec_clades[(b, m)].add(r["clade"])
        else:
            b = bag_of(r["clade"], bags)
            if b:
                real[b].append({**r, "_s": s})

    rows = []
    for bag in sorted(bags):
        L, S = lanes[bag], strata.get(bag, {})
        blks = sorted(real.get(bag, []), key=lambda r: -r["_s"])
        top = blks[0] if blks else None
        row = {
            "bag": bag, "lane": L["lane"], "lane_group": L["lane_group"],
            "top_domain": S.get("top_domain", ""),
            "top_domain_frac_of_reps": S.get("top_domain_frac_of_reps", ""),
            "all_domains": S.get("all_domains", ""),
            "n_distinct_anchors": L["n_distinct"], "family_size": L["family_size"],
            "step4_scored": bag in scored,
            "n_blocks": len(blks) if bag in scored else "",
            "max_composite_score": round(top["_s"], 4) if top else "",
            "top_block_id": top["block_id"] if top else "",
            "top_member_fraction": top.get("member_fraction", "") if top else "",
            "score_floor_used": floor if floor is not None else "",
            "n_blocks_above_floor": (
                sum(1 for b in blks if b["_s"] > floor)
                if floor is not None and bag in scored else ""),
        }
        for m in MODES:
            nb, ds = len(built[(bag, m)]), dec[(bag, m)]
            row[f"{m}__n_decoy_clades"] = nb
            row[f"{m}__n_decoy_blocks"] = len(ds)
            row[f"{m}__decoy_max"] = round(max(ds), 4) if ds else ""
            # available / clean / contaminated, never a bare True: a null that
            # built no decoys and a null whose decoys all stayed low are
            # opposite findings and must not share a cell value.
            row[f"{m}__status"] = (
                "not_available" if nb == 0 else
                "no_decoy_blocks" if not ds else
                "" if floor is None else
                "clean" if max(ds) <= floor else "contaminated")
        strictly_clean = all(
            row[f"{m}__status"] in ("clean", "no_decoy_blocks", "not_available")
            for m in STRICT)
        any_strict = any(row[f"{m}__n_decoy_clades"] for m in STRICT)
        row["strict_nulls_clean"] = strictly_clean and any_strict
        row["n_strict_decoy_clades"] = sum(
            row[f"{m}__n_decoy_clades"] for m in STRICT)
        row["floor_backed_by_n_decoy_blocks"] = floor_n
        rows.append(row)

    base = ["bag", "lane", "lane_group", "top_domain",
            "top_domain_frac_of_reps", "n_distinct_anchors", "family_size",
            "step4_scored", "n_blocks", "max_composite_score",
            "top_member_fraction", "score_floor_used", "n_blocks_above_floor",
            "strict_nulls_clean", "n_strict_decoy_clades",
            "floor_backed_by_n_decoy_blocks"]
    nullc = [f"{m}__{k}" for m in MODES
             for k in ("n_decoy_clades", "n_decoy_blocks", "decoy_max",
                       "status")]
    write(args.out / "all_bags_annotated.tsv",
          base + nullc + ["all_domains", "top_block_id"],
          sorted(rows, key=lambda r: (r["lane"], r["bag"])))

    # --- positive controls: recovery, not ranking -------------------------
    pos = [r for r in rows if r["lane"] == "known_rna_mge"]
    write(args.out / "positive_control_recovery.tsv", base + nullc,
          sorted(pos, key=lambda r: (r["lane_group"],
                                     -(r["max_composite_score"] or 0)
                                     if isinstance(r["max_composite_score"],
                                                   float) else 0)))
    recov = defaultdict(lambda: [0, 0, 0])
    for r in pos:
        g = recov[r["lane_group"]]
        g[0] += 1
        if r["step4_scored"] and r["n_blocks"]:
            g[1] += 1
        if r["strict_nulls_clean"] and isinstance(
                r["n_blocks_above_floor"], int) and r["n_blocks_above_floor"]:
            g[2] += 1

    # --- discovery: the result -------------------------------------------
    disc = [r for r in rows if r["lane"] in ("discovery", "unknown")
            and r["step4_scored"] and r["n_blocks"]]
    disc.sort(key=lambda r: (
        -int(bool(r["strict_nulls_clean"])),                       # 2
        -int(r[f"same_family_permutation__status"] == "clean"),    # 3
        -(r["max_composite_score"] if isinstance(
            r["max_composite_score"], float) else 0.0),            # 4
        -float(r["top_member_fraction"] or 0),                     # 5
        -NOVELTY.get(r["lane"], 0),                                # 6
        r["bag"]))
    write(args.out / "discovery_candidates.tsv", base + nullc, disc)

    # --- host-like: background pressure -----------------------------------
    host = [r for r in rows if r["lane"] == "host_like"]
    write(args.out / "host_like_background.tsv", base + nullc,
          sorted(host, key=lambda r: -(r["max_composite_score"]
                                       if isinstance(r["max_composite_score"],
                                                     float) else 0.0)))

    def rate(rs):
        run = [r for r in rs if r["step4_scored"]]
        hit = [r for r in run if r["n_blocks"]]
        return len(run), len(hit), (round(len(hit) / len(run), 3) if run else "")

    print(f"score floor (no observed decoy support) = {floor}  "
          f"backed by {floor_n} strict-null decoy blocks")
    print(f"legacy genomic floor was {fl.get('legacy_genomic_floor')}, "
          f"not used here\n")
    print(f"{'lane':16s} {'run':>5s} {'with block':>11s} {'rate':>6s}")
    for lane in ("known_rna_mge", "discovery", "unknown", "host_like"):
        n, h, p = rate([r for r in rows if r["lane"] == lane])
        print(f"{lane:16s} {n:5d} {h:11d} {str(p):>6s}")
    print("\n== positive control recovery (spike-in lane) ==")
    print(f"{'system':22s} {'bags':>5s} {'with block':>11s} "
          f"{'clean & above floor':>20s}")
    for g, (n, h, c) in sorted(recov.items()):
        print(f"{g:22s} {n:5d} {h:11d} {c:20d}")
    print(f"\n== top discovery candidates (known RNA-MGE excluded) ==")
    for r in disc[:20]:
        print(f"  {r['bag']}  {r['lane']:10s} {str(r['max_composite_score']):>7s} "
              f"clean={r['strict_nulls_clean']!s:5s} "
              f"n={r['n_strict_decoy_clades']:>3d}  {r['top_domain']}")
    print(f"\n-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
