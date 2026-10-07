#!/usr/bin/env python
"""Locus-level triage for the step-4 candidate blocks.

    triage_report.py --candidates <dir>/step4_candidates.tsv \\
                     --blocks <k3_dir> --windows <windows_dir> \\
                     --fold <fold_dir> --padding <fasta_dir>/block_padding.tsv \\
                     --out <dir> [--exclude-benchmark-only]

One row per block (`triage_blocks.tsv`) and one per locus group
(`triage_loci.tsv`). The point of the grouping is that a block is not an
element: a structured RNA interrupted by its own ORF, or by a prodigal call that
lands inside it, is reported by step 4 as two or more blocks. Reading them
separately understates the locus and overstates the block count.

Blocks of one clade are merged into a `locus_group_id` when either

  * both lie within `--merge-window` bp of the anchor ORF, or
  * they abut opposite sides of it (one ending at the start codon, one
    beginning at the stop),

by union-find, so a chain of blocks around one anchor becomes one group. The
rule is anchor-geometric and clade-local; it never looks at sequence, and never
at a known RNA.

Blind: coordinates, ORF geometry and each block's own fold. Benchmark columns
are annotation only -- they are joined on after every other column is computed,
and nothing here reads them back.
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_ncrna.compute import require_slurm

_d = __import__("importlib.util", fromlist=["util"])
_spec = _d.spec_from_file_location("_db", Path(__file__).resolve().parent / "discover_blocks.py")
_db = _d.module_from_spec(_spec)
_spec.loader.exec_module(_db)

# R-scape's own estimate of how many pairs it *could* have seen covary. Below 1
# the test never had a chance; below 3 a single pair carries the verdict.
POWER_NONE, POWER_LOW = 1.0, 3.0


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        cols = next(fh).rstrip("\n").split("\t")
        return [dict(zip(cols, line.rstrip("\n").split("\t"))) for line in fh]


def write_tsv(path: Path, rows: list[dict], cols: list[str]) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def q(v: list[float], f: float):
    s = sorted(v)
    return s[min(len(s) - 1, int(f * len(s)))]


def num(s, d=None):
    try:
        return float(s)
    except (TypeError, ValueError):
        return d


class UF:
    def __init__(self):
        self.p: dict[str, str] = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def side_of(rows: list[dict]) -> str:
    """Which side of the anchor ORF the members sit on.

    `split` is a real answer, not a fallback: a block whose members land
    upstream in some genomes and downstream in others is a different object
    from one that is consistently placed, and merging it into a locus on
    median coordinates alone would hide that.

    Each member is tested against *its own* anchor length. Using the clade
    median instead manufactures split calls out of nothing: a block that starts
    exactly at the stop codon in every genome is downstream in every genome, but
    against a median anchor length every member with a slightly longer anchor
    scores as overlapping. That is what turned a 120/120-downstream IS110 block
    into a `split` call.
    """
    up = sum(1 for r in rows if r["end"] <= 0)
    dn = sum(1 for r in rows if r["start"] >= r["anchor_len"])
    ov = len(rows) - up - dn
    n = max(1, len(rows))
    for name, c in (("upstream", up), ("downstream", dn), ("overlapping", ov)):
        if c / n >= 0.9:
            return name
    return "split"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True, type=Path)
    ap.add_argument("--blocks", required=True, type=Path)
    ap.add_argument("--windows", required=True, type=Path)
    ap.add_argument("--fold", type=Path)
    ap.add_argument("--padding", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--merge-window", type=int, default=500)
    ap.add_argument("--exclude-benchmark-only", action="store_true")
    args = ap.parse_args()
    require_slurm("triage report")
    args.out.mkdir(parents=True, exist_ok=True)

    cands = read_tsv(args.candidates)
    if args.exclude_benchmark_only:
        cands = [c for c in cands if c["benchmark_only"] != "True"]
    by_id = {c["block_id"]: c for c in cands}
    print(f"{len(cands)} candidate blocks")

    blocks = {b["block_id"]: b for b in read_tsv(args.blocks / "blocks.tsv")}
    members = defaultdict(list)
    for m in read_tsv(args.blocks / "block_members.tsv"):
        if m["block_id"] in by_id:
            members[m["block_id"]].append(m)

    fold = {f["block_id"]: f for f in read_tsv(args.fold / "fold_summary.tsv")} \
        if args.fold else {}
    padd = {p["block_id"]: p for p in read_tsv(args.padding)} if args.padding else {}

    # Anchor lengths come from the window headers; no prodigal run is needed
    # here because nothing in this report depends on the coding mask.
    anchor_len: dict[str, dict[str, int]] = {}
    for clade in {c["clade_id"] for c in cands}:
        anchor_len[clade] = {r["anchor_id"]: r["anchor_len"]
                             for r in _db.load_clade(args.windows / f"{clade}.fasta")}

    rows = []
    for c in cands:
        bid, clade = c["block_id"], c["clade_id"]
        al = anchor_len.get(clade, {})
        ms = [{"start": int(m["start"]), "end": int(m["end"]),
               "anchor_len": al.get(m["anchor_id"], 0)}
              for m in members.get(bid, []) if m["anchor_id"] in al]
        if not ms:
            print(f"  {bid}: no members resolved", file=sys.stderr)
            continue
        S = [m["start"] for m in ms]
        E = [m["end"] for m in ms]
        a_len = int(statistics.median([m["anchor_len"] for m in ms]))
        f = fold.get(bid, {})
        p = padd.get(bid, {})
        exp = num(f.get("expected_covarying"))
        tp, fp = num(f.get("covarying_pairs")), num(f.get("false_positive_pairs"))
        # R-scape omits its power section entirely when an alignment has no
        # substitutions at all (it reports avgid 100). That is the extreme of
        # "none", not a missing measurement, and must not read as blank.
        avgid = num(f.get("rscape_avgid"))
        power = ("none_no_variation" if exp is None and avgid is not None and avgid >= 99.0
                 else "unknown" if exp is None
                 else "none" if exp < POWER_NONE
                 else "low" if exp < POWER_LOW else "adequate")
        rows.append({
            "family": c["family"], "clade_id": clade, "block_id": bid,
            "score": c["score"],
            "covariation_verdict": f.get("verdict", "not_folded"),
            "rscape_power": power,
            "rscape_expected_covarying": f.get("expected_covarying", ""),
            "rscape_bpairs_tested": f.get("bpairs_tested", ""),
            "TP_covarying_pairs": f.get("covarying_pairs", ""),
            "FP_pairs": f.get("false_positive_pairs", ""),
            "PPV": f.get("ppv", ""),
            "member_fraction": c["member_fraction"], "n_members": len(ms),
            "relative_start_median": int(statistics.median(S)),
            "relative_start_q25": int(q(S, .25)), "relative_start_q75": int(q(S, .75)),
            "relative_end_median": int(statistics.median(E)),
            "relative_end_q25": int(q(E, .25)), "relative_end_q75": int(q(E, .75)),
            "side_of_anchor": side_of(ms),
            "distance_to_anchor": int(blocks[bid]["dist_to_anchor"]),
            "anchor_len": a_len,
            "padded_coding_frac": p.get("padded_coding_frac", ""),
            "median_pairwise_identity": f.get("mean_pairwise_identity", ""),
            # Annotation only, joined after everything above is fixed.
            "benchmark_overlaps_known_ncrna": c.get("benchmark_overlaps_known_ncrna", ""),
            "benchmark_max_overlap_bp": c.get("benchmark_max_overlap_bp", ""),
        })
    by_block = {r["block_id"]: r for r in rows}

    # --- neighbours and locus grouping ------------------------------------- #
    uf, nf = UF(), UF()
    for r in rows:
        uf.find(r["block_id"])
        nf.find(r["block_id"])
    per_clade = defaultdict(list)
    for r in rows:
        per_clade[r["clade_id"]].append(r)
    for clade, rs in per_clade.items():
        for i, a in enumerate(rs):
            nb = []
            for b in rs:
                if b is a:
                    continue
                # Gap between the two blocks' median intervals, 0 if they touch.
                gap = max(a["relative_start_median"] - b["relative_end_median"],
                          b["relative_start_median"] - a["relative_end_median"], 0)
                if gap <= args.merge_window:
                    nb.append(b["block_id"].rsplit("__", 1)[-1])
                    nf.union(a["block_id"], b["block_id"])
                near_anchor = (abs(a["distance_to_anchor"]) <= args.merge_window
                               and abs(b["distance_to_anchor"]) <= args.merge_window)
                opposite = ({a["side_of_anchor"], b["side_of_anchor"]}
                            == {"upstream", "downstream"}
                            and abs(a["distance_to_anchor"]) <= args.merge_window
                            and abs(b["distance_to_anchor"]) <= args.merge_window)
                if near_anchor or opposite:
                    uf.union(a["block_id"], b["block_id"])
            a["same_clade_neighbor_blocks_500bp"] = ",".join(sorted(nb)) or "-"

    # Stable, readable group ids: family + clade ordinal + group ordinal.
    groups = defaultdict(list)
    for r in rows:
        groups[uf.find(r["block_id"])].append(r)
    clade_ord, gid_of = {}, {}
    for root, rs in sorted(groups.items(),
                           key=lambda kv: (kv[1][0]["family"], kv[1][0]["clade_id"],
                                           -max(float(x["score"]) for x in kv[1]))):
        cl = rs[0]["clade_id"]
        clade_ord.setdefault(cl, len(clade_ord) + 1)
        n = sum(1 for g in gid_of.values() if g.startswith(
            f"{rs[0]['family']}_C{clade_ord[cl]:02d}_"))
        gid_of[root] = f"{rs[0]['family']}_C{clade_ord[cl]:02d}_L{n + 1}"
    for r in rows:
        r["locus_group_id"] = gid_of[uf.find(r["block_id"])]

    # `nearby_block_group_id` is block-to-block proximity within a clade, and
    # nothing more. It is deliberately *not* folded into locus_group_id: the
    # question this pipeline asks is whether a novel protein ORF has an ncRNA
    # beside it, and two blocks sitting >1 kb from the anchor may be one distal
    # RNA or may be the UTR and terminator of some other gene. Merging them into
    # a protein-associated locus would blur that distinction; flagging them
    # keeps the adjacency visible without asserting it means anything.
    nsize: dict[str, int] = defaultdict(int)
    for r in rows:
        nsize[nf.find(r["block_id"])] += 1
    nid, seen = {}, 0
    for r in sorted(rows, key=lambda r: (r["clade_id"], r["relative_start_median"])):
        root = nf.find(r["block_id"])
        if nsize[root] < 2:
            r["nearby_block_group_id"] = "-"
            continue
        if root not in nid:
            seen += 1
            nid[root] = f"{r['family']}_N{seen:02d}"
        r["nearby_block_group_id"] = nid[root]

    bcols = ["locus_group_id", "family", "clade_id", "block_id", "score",
             "covariation_verdict", "rscape_power", "rscape_expected_covarying",
             "rscape_bpairs_tested", "TP_covarying_pairs", "FP_pairs", "PPV",
             "member_fraction", "n_members", "relative_start_median",
             "relative_start_q25", "relative_start_q75", "relative_end_median",
             "relative_end_q25", "relative_end_q75", "side_of_anchor",
             "distance_to_anchor", "anchor_len", "same_clade_neighbor_blocks_500bp",
             "nearby_block_group_id",
             "padded_coding_frac", "median_pairwise_identity",
             "benchmark_overlaps_known_ncrna", "benchmark_max_overlap_bp"]
    rows.sort(key=lambda r: (r["locus_group_id"], -float(r["score"])))
    write_tsv(args.out / "triage_blocks.tsv", rows, bcols)

    # --- one row per locus -------------------------------------------------- #
    lrows = []
    for gid in sorted({r["locus_group_id"] for r in rows}):
        rs = sorted((r for r in rows if r["locus_group_id"] == gid),
                    key=lambda r: -float(r["score"]))
        sides = {r["side_of_anchor"] for r in rs}
        arch = ("flanks_both_sides" if {"upstream", "downstream"} <= sides else
                "split" if "split" in sides else
                "+".join(sorted(sides)))
        tps = [num(r["TP_covarying_pairs"], 0) for r in rs]
        ids = [num(r["median_pairwise_identity"]) for r in rs]
        ids = [i for i in ids if i is not None]
        ver = ["covariation_supported", "covariation_borderline",
               "tested_no_covariation", "underpowered", "covariation_test_failed",
               "not_folded"]
        best = next((v for v in ver if v in {r["covariation_verdict"] for r in rs}),
                    "")
        lrows.append({
            "locus_group_id": gid, "family": rs[0]["family"],
            "clade_id": rs[0]["clade_id"], "n_blocks": len(rs),
            "block_ids": ",".join(r["block_id"].rsplit("__", 1)[-1] for r in rs),
            "max_score": rs[0]["score"],
            "best_covariation_verdict": best,
            "n_blocks_covariation_supported":
                sum(1 for r in rs if r["covariation_verdict"] == "covariation_supported"),
            "n_blocks_covariation_borderline":
                sum(1 for r in rs if r["covariation_verdict"] == "covariation_borderline"),
            "total_TP_covarying_pairs": int(sum(tps)),
            "total_FP_pairs": int(sum(num(r["FP_pairs"], 0) for r in rs)),
            "architecture": arch,
            "locus_start_median": min(r["relative_start_median"] for r in rs),
            "locus_end_median": max(r["relative_end_median"] for r in rs),
            "locus_span_bp": (max(r["relative_end_median"] for r in rs)
                              - min(r["relative_start_median"] for r in rs)),
            "block_bp_in_locus": sum(r["relative_end_median"] - r["relative_start_median"]
                                     for r in rs),
            "anchor_len": rs[0]["anchor_len"],
            "max_n_members": max(r["n_members"] for r in rs),
            "max_member_fraction": max(float(r["member_fraction"]) for r in rs),
            "median_pairwise_identity": round(statistics.median(ids), 4) if ids else "",
            "nearby_block_group_ids": ",".join(sorted(
                {r["nearby_block_group_id"] for r in rs} - {"-"})) or "-",
            "benchmark_any_known_overlap":
                any(r["benchmark_overlaps_known_ncrna"] == "True" for r in rs),
        })
    lcols = ["locus_group_id", "family", "clade_id", "n_blocks", "block_ids",
             "max_score", "best_covariation_verdict", "n_blocks_covariation_supported",
             "n_blocks_covariation_borderline",
             "total_TP_covarying_pairs", "total_FP_pairs", "architecture",
             "locus_start_median", "locus_end_median", "locus_span_bp",
             "block_bp_in_locus", "anchor_len", "max_n_members",
             "max_member_fraction", "median_pairwise_identity",
             "nearby_block_group_ids",
             "benchmark_any_known_overlap"]
    lrows.sort(key=lambda r: (-float(r["max_score"])))
    write_tsv(args.out / "triage_loci.tsv", lrows, lcols)

    summary = {
        "n_blocks": len(rows), "n_loci": len(lrows),
        "merge_window_bp": args.merge_window,
        "loci_with_multiple_blocks": sum(1 for r in lrows if r["n_blocks"] > 1),
        "architecture": {a: sum(1 for r in lrows if r["architecture"] == a)
                         for a in sorted({r["architecture"] for r in lrows})},
        "verdicts": {v: sum(1 for r in rows if r["covariation_verdict"] == v)
                     for v in sorted({r["covariation_verdict"] for r in rows})},
        "power": {p: sum(1 for r in rows if r["rscape_power"] == p)
                  for p in sorted({r["rscape_power"] for r in rows})},
        "loci_with_support": sum(1 for r in lrows
                                 if r["n_blocks_covariation_supported"] > 0),
        "loci_borderline_only": sum(
            1 for r in lrows if r["n_blocks_covariation_supported"] == 0
            and r["best_covariation_verdict"] == "covariation_borderline"),
        "nearby_block_groups": len({r["nearby_block_group_id"] for r in rows} - {"-"}),
    }
    (args.out / "triage.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    w = ["locus_group_id", "architecture", "max_score", "n_blocks",
         "total_TP_covarying_pairs", "best_covariation_verdict"]
    print("\n" + "  ".join(w))
    for r in lrows:
        print("  ".join(str(r[c]).ljust(len(c)) for c in w))
    print(f"\n-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
