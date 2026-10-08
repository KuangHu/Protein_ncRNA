#!/usr/bin/env python
"""UNBLIND interpretation layer for a single step-4 block.

    interpret_locus.py --block <BLOCK_ID> --run final_26739540 \\
                       --scratch /global/scratch/users/kh36969/protein_ncrna

Everything upstream of this file is blind: no curated RNA entered anchor
selection, clade decomposition, window extraction, block calling, scoring or
folding. This script is where that stops. It asks what a block *is* -- does it
overlap a RetronDB ncRNA, do defence or mobile-element markers sit beside its
anchors, do its covarying pairs land on real helices -- and to do that it reads
annotation.

Keeping it in a separate file with `UNBLIND` in the name is the point. The
discovery score was computed before any of this was consulted, and nothing here
can feed back into it: the script reads frozen outputs and writes a report. If a
future change makes a column from here reachable from a scoring path,
`tests/test_leakage.py` should fail, and that is the intended behaviour.

Sections written to `<out>/<short>.interpretation.md` and `.json`:
  1. clade composition -- anchors, distinct proteins, species, identity
  2. block geometry -- from the blind locus table, restated for context
  3. known-RNA overlap -- RetronDB, per block and per clade
  4. neighbourhood markers -- cascade / HNH / TOPRIM / TIR beside the anchors
  5. covariation -- covarying pairs grouped into the helices they support
"""

from __future__ import annotations

import argparse
import json
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

# Distance within which a covarying pair's two positions are called the same
# helix. Pairs stack, so consecutive pairs of a helix differ by ~1 on each side.
HELIX_GAP = 4


def read_tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with open(path) as fh:
        cols = next(fh).rstrip("\n").split("\t")
        return [dict(zip(cols, line.rstrip("\n").split("\t"))) for line in fh]


def covarying_pairs(cov: Path) -> list[dict]:
    out = []
    if not cov.exists():
        return out
    for line in open(cov):
        if line.startswith(("#", "\n")):
            continue
        f = line.split()
        if len(f) < 8 or not f[1].isdigit():
            continue
        out.append({"i": int(f[1]), "j": int(f[2]), "score": float(f[3]),
                    "evalue": float(f[4]), "substitutions": int(f[6]),
                    "power": float(f[7])})
    return sorted(out, key=lambda p: p["i"])


def group_helices(pairs: list[dict]) -> list[dict]:
    """Collapse stacked covarying pairs into the helices they support.

    14 covarying pairs spread over 7 helices is a different claim from 14 pairs
    in one hairpin: the first says the fold is right in several places, the
    second says it is right in one. R-scape reports pairs, so the grouping has
    to happen here.
    """
    helices: list[dict] = []
    for p in pairs:
        for h in helices:
            if (abs(p["i"] - h["i_last"]) <= HELIX_GAP
                    and abs(p["j"] - h["j_last"]) <= HELIX_GAP):
                h["n"] += 1
                h["i_last"], h["j_last"] = p["i"], p["j"]
                h["i_end"], h["j_end"] = p["i"], p["j"]
                h["max_power"] = max(h["max_power"], p["power"])
                h["min_evalue"] = min(h["min_evalue"], p["evalue"])
                break
        else:
            helices.append({"i_start": p["i"], "j_start": p["j"],
                            "i_end": p["i"], "j_end": p["j"],
                            "i_last": p["i"], "j_last": p["j"], "n": 1,
                            "max_power": p["power"], "min_evalue": p["evalue"]})
    for h in helices:
        h.pop("i_last"), h.pop("j_last")
        h["loop_span"] = h["j_start"] - h["i_end"]
    return helices


def structure_string(path: Path) -> tuple[str, str]:
    """(comment header, dot-bracket consensus) from a fold_blocks structure.txt."""
    if not path.exists():
        return "", ""
    lines = [l.rstrip("\n") for l in open(path)]
    head = [l for l in lines if l.startswith("#")]
    ss = "".join(l for l in lines if l and not l.startswith("#"))
    return "\n".join(head), ss


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--block", required=True)
    ap.add_argument("--run", default="final_26739540")
    ap.add_argument("--scratch", type=Path,
                    default=Path("/global/scratch/users/kh36969/protein_ncrna"))
    ap.add_argument("--neighborhood", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    S, run = args.scratch, args.run
    cand = S / f"step4_candidates_{run}"
    win = S / f"windows_{run}"
    bench = S / f"known_retron_benchmark_{run}"
    nbr = args.neighborhood or (S / "nbr_26710278" / "neighborhood.tsv")
    out = args.out or (cand / "interpretation")
    out.mkdir(parents=True, exist_ok=True)

    bid = args.block
    clade = bid.rsplit("__", 1)[0]
    fam = bid.split("__", 1)[0]
    short = f"{fam}_{bid.rsplit('__', 1)[-1]}"
    m = re.search(r"GC[AF]_(\d+)", bid)
    if m:
        short = f"{fam}_{m.group(1)[:10]}_{bid.rsplit('__', 1)[-1]}"

    R: dict = {"block_id": bid, "clade_id": clade, "family": fam, "run": run}

    # --- 1. clade composition -------------------------------------------
    cap = {r["clade_id"]: r for r in read_tsv(win / "cap_diagnostics.tsv")}
    R["clade"] = {k: v for k, v in (cap.get(clade) or {}).items()
                  if k != "clade_id"}
    for row in read_tsv(S / f"clades_{run}" / "selected_clades.tsv"):
        if row.get("cluster_id", "").replace("|", "_").replace(":", "_") in clade \
                or row.get("cluster_id") in clade:
            R["clade"].update({k: row[k] for k in
                               ("n_anchors", "median_identity", "n_species",
                                "median_aa_len", "n_complete_windows")
                               if k in row})
            break

    # --- 2. block geometry (restated from the blind triage) --------------
    for row in read_tsv(cand / "triage" / "triage_blocks.tsv"):
        if row["block_id"] == bid:
            R["geometry"] = {k: row[k] for k in (
                "score", "side_of_anchor", "relative_start_median",
                "relative_start_q25", "relative_start_q75",
                "relative_end_median", "relative_end_q25", "relative_end_q75",
                "anchor_len", "n_members", "member_fraction",
                "padded_coding_frac", "median_pairwise_identity",
                "covariation_verdict", "rscape_power",
                "same_clade_neighbor_blocks_500bp") if k in row}
            break

    # --- 3. known-RNA overlap (RetronDB) ---------------------------------
    ov = {}
    for row in read_tsv(bench / "benchmark_block_overlap.tsv"):
        if row["block_id"] == bid:
            ov = {"overlaps_known_retron": row["overlaps_known_retron"],
                  "members_overlapping_known": row["members_overlapping_known"],
                  "max_overlap_bp": row["max_overlap_bp"],
                  "rank_in_clade": row["rank_in_clade"]}
            break
    for row in read_tsv(bench / "benchmark_clade_overlap.tsv"):
        if row.get("cluster_id", "") and row["cluster_id"].replace("|", "_") in \
                clade.replace("|", "_"):
            ov["clade_known_retron_anchors"] = row.get("n_known_retron_anchors")
            ov["clade_known_frac"] = row.get("known_frac")
            ov["clade_retron_subtypes"] = row.get("retron_subtypes")
            ov["clade_top_retron_names"] = row.get("top_retron_names")
            break
    R["known_rna_overlap"] = ov or {"overlaps_known_retron": "False",
                                    "note": "block not present in benchmark table"}

    # --- 4. neighbourhood markers ----------------------------------------
    members = [r for r in read_tsv(S / f"blocks_{run}" / "block_members.tsv")
               if r.get("block_id") == bid]
    anchors = {r["anchor_id"] for r in members if "anchor_id" in r}
    marks: dict[str, list[int]] = defaultdict(list)
    if nbr.exists() and anchors:
        for row in read_tsv(nbr):
            if row["anchor_id"] in anchors:
                marks[row["neighbor_group"]].append(int(row["distance_bp"]))
    R["neighborhood"] = {
        "n_block_members": len(members),
        "n_anchors_checked": len(anchors),
        "markers": {g: {"n_anchors": len(set(d)), "n_hits": len(d),
                        "median_distance_bp": int(statistics.median(d))}
                    for g, d in sorted(marks.items())},
        "interpretation_note": (
            "Cascade/Cas1/Cas5/Cas6 beside an anchor supports a CRISPR locus; "
            "HNH/TOPRIM/TIR/HEPN beside an RT supports a defence RT without "
            "naming it a retron, which RT homology alone cannot establish. An "
            "empty marker set is not evidence of absence -- only these groups "
            "were scanned."),
    }

    # --- 5. covariation ---------------------------------------------------
    fold = cand / "fold_all080" / bid
    head, ss = structure_string(fold / "structure.txt")
    pairs = covarying_pairs(fold / "rscape" / "aln_1.cov")
    hel = group_helices(pairs)
    R["covariation"] = {
        "n_covarying_pairs": len(pairs),
        "n_helices_supported": len(hel),
        "alignment_len": len(ss),
        "consensus_header": head,
        "helices": hel,
        "pairs": pairs,
    }

    # --- write -------------------------------------------------------------
    (out / f"{short}.interpretation.json").write_text(json.dumps(R, indent=2))

    L = [f"# Interpretation: `{bid}`", "",
         "**UNBLIND.** Everything above step 5 was computed without reading any",
         "curated RNA. This file reads annotation and cannot feed back into the",
         "score.", "",
         "## 1. Clade composition", ""]
    for k, v in R["clade"].items():
        L.append(f"- `{k}`: {v}")
    L += ["", "## 2. Block geometry (restated from the blind triage)", ""]
    for k, v in R.get("geometry", {}).items():
        L.append(f"- `{k}`: {v}")
    L += ["", "## 3. Known-RNA overlap (RetronDB)", ""]
    for k, v in R["known_rna_overlap"].items():
        L.append(f"- `{k}`: {v}")
    L += ["", "## 4. Neighbourhood markers", "",
          f"{R['neighborhood']['n_anchors_checked']} anchors checked."]
    if R["neighborhood"]["markers"]:
        L += ["", "| marker | anchors | median distance bp |", "|---|---:|---:|"]
        for g, d in R["neighborhood"]["markers"].items():
            L.append(f"| {g} | {d['n_anchors']} | {d['median_distance_bp']} |")
    else:
        L.append("\nNo scanned marker within the scan window of any member anchor.")
    L += ["", f"_{R['neighborhood']['interpretation_note']}_", "",
          "## 5. Covariation", "",
          f"{len(pairs)} covarying pairs over {len(hel)} helices, "
          f"alignment length {len(ss)}.", ""]
    if hel:
        L += ["| helix | i | j | pairs | loop span | max power | min E |",
              "|---:|---:|---:|---:|---:|---:|---:|"]
        for n, h in enumerate(hel, 1):
            L.append(f"| {n} | {h['i_start']}–{h['i_end']} | "
                     f"{h['j_end']}–{h['j_start']} | {h['n']} | "
                     f"{h['loop_span']} | {h['max_power']:.2f} | "
                     f"{h['min_evalue']:.2g} |")
    (out / f"{short}.interpretation.md").write_text("\n".join(L) + "\n")

    print(json.dumps({k: v for k, v in R.items() if k != "covariation"}, indent=2))
    print(f"\ncovariation: {len(pairs)} pairs over {len(hel)} helices")
    for n, h in enumerate(hel, 1):
        print(f"  helix {n}: {h['i_start']}-{h['i_end']} x {h['j_end']}-{h['j_start']}"
              f"  pairs={h['n']} loop={h['loop_span']} power={h['max_power']:.2f}")
    print(f"-> {out}/{short}.interpretation.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
