#!/usr/bin/env python
"""A rank must not turn an untested candidate into a refuted one.

The homology-pool run produced three outcomes that a pass/fail column flattens
into the same word: a null that was built and lost, a null that was built and
produced nothing at all, and a null that could not be built. The ranking exists
to keep them apart, and these tests pin the two ways it could stop doing that.

The failure that matters most is CDS07978's: the highest-recurrence block in the
corpus, in a bag with no homolog anywhere and no Pfam family, so no hard null is
constructible. If the ranking treats "untested" as "clean" it leads the entire
table on evidence it never had; if it treats "untested" as "failed" it buries
the single most interesting thing the run found. It must do neither, and the
information gap is what carries the distinction.
"""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable
ENV = {"PN_ALLOW_LOGIN_NODE": "1", "PATH": "/usr/bin:/bin"}

CAND_COLS = ["bag", "lane", "n_distinct_anchors", "n_blocks",
             "max_composite_score", "n_blocks_above_floor"]
NULL_COLS = ["bag", "null_type", "n_decoy_clades", "n_decoy_blocks",
             "decoy_max_score"]
POOL_COLS = ["target", "null_strength_source", "n_source_bags",
             "pool_distinct_proteins"]


def _tsv(p: Path, cols, rows):
    with open(p, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def _run(tmp: Path, cand, null, pool, floor=0.7077):
    _tsv(tmp / "cand.tsv", CAND_COLS, cand)
    _tsv(tmp / "null.tsv", NULL_COLS, null)
    _tsv(tmp / "pool.tsv", POOL_COLS, pool)
    (tmp / "targets.txt").write_text("\n".join(c["bag"] for c in cand) + "\n")
    r = subprocess.run(
        [PY, "-I", str(ROOT / "scripts" / "rank_candidates.py"),
         "--candidates", str(tmp / "cand.tsv"), "--null", str(tmp / "null.tsv"),
         "--pool", str(tmp / "pool.tsv"), "--targets", str(tmp / "targets.txt"),
         "--floor", str(floor), "--out", str(tmp / "out")],
        capture_output=True, text=True, env=ENV)
    assert r.returncode == 0, r.stderr
    return {x["bag"]: x for x in csv.DictReader(
        open(tmp / "out" / "ranked_candidates.tsv"), delimiter="\t")}


def test_an_untested_block_outranks_nothing_on_strength_it_never_showed(tmp_path):
    """Equal recurrence, equal breadth; one faced a working null and won.

    The tested one must rank above the untested one. Ranking on recurrence alone
    would tie them, which is the whole error: surviving a null that demonstrably
    produces blocks is evidence, and never having been asked is not.
    """
    cand = [
        {"bag": "TESTED", "lane": "unknown", "n_distinct_anchors": 400,
         "n_blocks": 3, "max_composite_score": 0.95, "n_blocks_above_floor": 2},
        {"bag": "UNTESTED", "lane": "unknown", "n_distinct_anchors": 400,
         "n_blocks": 3, "max_composite_score": 0.95, "n_blocks_above_floor": 2},
    ]
    null = [
        {"bag": "TESTED", "null_type": "homology_pool_permutation",
         "n_decoy_clades": 6, "n_decoy_blocks": 2, "decoy_max_score": 0.60},
        {"bag": "UNTESTED", "null_type": "homology_pool_permutation",
         "n_decoy_clades": 0, "n_decoy_blocks": 0, "decoy_max_score": ""},
    ]
    pool = [{"target": "TESTED", "null_strength_source": "strong"},
            {"target": "UNTESTED", "null_strength_source": "weak"}]
    out = _run(tmp_path, cand, null, pool)

    assert int(out["TESTED"]["rank"]) < int(out["UNTESTED"]["rank"])
    assert out["TESTED"]["tier"] == "A"
    assert out["UNTESTED"]["tier"] == "D"
    # Identical signal: the ranking must separate them on the test, not by
    # quietly discounting the untested block's quality.
    assert out["TESTED"]["signal"] == out["UNTESTED"]["signal"]
    # ... and the untested one keeps a real score. Zero would read as refuted.
    assert float(out["UNTESTED"]["evidence_score"]) > 0
    # The whole gap is attributable to the missing test.
    assert float(out["UNTESTED"]["information_gap"]) > float(
        out["TESTED"]["information_gap"])
    assert float(out["TESTED"]["information_gap"]) == 0.0


def test_a_near_miss_is_demoted_below_a_quiet_null(tmp_path):
    """A decoy that came within 0.06 is the worst news in the table.

    Nominally both bags are `clean` -- no decoy crossed the floor. But one had a
    null that produced a block almost as good as the real one, and that is
    evidence against; the other's null produced nothing. If the near miss does
    not fall below, a working null counts for less than an inert one.
    """
    cand = [
        {"bag": "NEARMISS", "lane": "unknown", "n_distinct_anchors": 300,
         "n_blocks": 3, "max_composite_score": 0.77, "n_blocks_above_floor": 2},
        {"bag": "QUIET", "lane": "unknown", "n_distinct_anchors": 300,
         "n_blocks": 3, "max_composite_score": 0.88, "n_blocks_above_floor": 2},
    ]
    null = [
        {"bag": "NEARMISS", "null_type": "homology_pool_permutation",
         "n_decoy_clades": 6, "n_decoy_blocks": 3, "decoy_max_score": 0.7077},
        {"bag": "QUIET", "null_type": "homology_pool_permutation",
         "n_decoy_clades": 3, "n_decoy_blocks": 0, "decoy_max_score": ""},
    ]
    pool = [{"target": "NEARMISS", "null_strength_source": "weak"},
            {"target": "QUIET", "null_strength_source": "strong"}]
    out = _run(tmp_path, cand, null, pool)

    assert out["NEARMISS"]["tier"] == "B"
    assert out["QUIET"]["tier"] == "C"
    assert int(out["QUIET"]["rank"]) < int(out["NEARMISS"]["rank"])
    # The margin is measured against the decoy, not the floor, precisely because
    # a decoy got there first.
    assert out["NEARMISS"]["challenge_basis"] == "decoy_max"
    assert out["QUIET"]["challenge_basis"] == "floor"
    assert float(out["NEARMISS"]["margin"]) < 0.10


def test_breadth_separates_equal_recurrence_over_unequal_lineages(tmp_path):
    """0.99 across 400 distinct anchors is not 0.99 across 25.

    composite_score is very nearly recurrence alone for these bags -- the
    geometric terms all saturate at 1.0 -- so without a breadth term the ranking
    has only one axis and a block seen in 25 near-identical lineages ties one
    seen in 400.
    """
    cand = [
        {"bag": "BROAD", "lane": "unknown", "n_distinct_anchors": 400,
         "n_blocks": 1, "max_composite_score": 0.90, "n_blocks_above_floor": 1},
        {"bag": "NARROW", "lane": "unknown", "n_distinct_anchors": 25,
         "n_blocks": 1, "max_composite_score": 0.90, "n_blocks_above_floor": 1},
    ]
    null = [{"bag": b, "null_type": "homology_pool_permutation",
             "n_decoy_clades": 3, "n_decoy_blocks": 0, "decoy_max_score": ""}
            for b in ("BROAD", "NARROW")]
    pool = [{"target": "BROAD"}, {"target": "NARROW"}]
    out = _run(tmp_path, cand, null, pool)

    assert float(out["BROAD"]["signal"]) > float(out["NARROW"]["signal"])
    assert int(out["BROAD"]["rank"]) < int(out["NARROW"]["rank"])
    # Narrow is discounted, not erased -- it is still a real recurrent block.
    assert float(out["NARROW"]["breadth_factor"]) >= 0.5


def test_only_the_named_null_mode_is_read(tmp_path):
    """within_bag_shuffle produced no block anywhere in this run.

    Its rows sit in the same table as the homology-pool rows. Reading them would
    report every bag as facing a null that never challenged anything, which is
    how `strict_nulls_clean = True` came to be set for a bag with no pool at all.
    """
    cand = [{"bag": "A", "lane": "unknown", "n_distinct_anchors": 400,
             "n_blocks": 1, "max_composite_score": 0.9,
             "n_blocks_above_floor": 1}]
    null = [
        {"bag": "A", "null_type": "within_bag_shuffle",
         "n_decoy_clades": 9, "n_decoy_blocks": 0, "decoy_max_score": ""},
        {"bag": "A", "null_type": "homology_pool_permutation",
         "n_decoy_clades": 0, "n_decoy_blocks": 0, "decoy_max_score": ""},
    ]
    out = _run(tmp_path, cand, null, [{"target": "A"}])
    assert out["A"]["tier"] == "D"
    assert out["A"]["n_decoy_clades"] == "0"


def test_an_orphan_is_not_the_same_as_having_no_curated_domain(tmp_path):
    """Four of the five domain-free bags have sub-threshold Pfam leads.

    Only CDS07978 has nothing even with the gathering cutoff removed. Collapsing
    the two would call five bags orphans when one is, and an orphan carrying a
    highly recurrent block is the most interesting object in the run.
    """
    dom = tmp_path / "t.domtbl"
    cols = ["X"] * 23
    def row(target, query, ie, sc, desc):
        f = list(cols)
        f[0], f[3], f[12], f[13], f[22] = target, query, str(ie), str(sc), desc
        return " ".join(f)
    dom.write_text("\n".join([
        "# comment",
        row("DUF1672", "HASLEAD|0", 0.0042, 16.8, "Protein of unknown function"),
        row("O_anti_polymase", "ORPHAN|0", 0.24, 11.5, "Putative O-antigen"),
        row("rve", "CURATED|0", 1e-30, 200.0, "Integrase core domain"),
    ]) + "\n")
    strict = tmp_path / "strict.tsv"
    _tsv(strict, ["bag", "all_domains"],
         [{"bag": "HASLEAD", "all_domains": ""},
          {"bag": "ORPHAN", "all_domains": ""},
          {"bag": "CURATED", "all_domains": "rve"}])
    (tmp_path / "tg.txt").write_text("HASLEAD\nORPHAN\nCURATED\n")
    r = subprocess.run(
        [PY, "-I", str(ROOT / "scripts" / "pfam_sensitive_report.py"),
         "--domtbl", str(dom), "--strict", str(strict),
         "--targets", str(tmp_path / "tg.txt"),
         "--out", str(tmp_path / "ps")],
        capture_output=True, text=True, env=ENV)
    assert r.returncode == 0, r.stderr
    got = {x["bag"]: x["sensitive_status"] for x in csv.DictReader(
        open(tmp_path / "ps" / "pfam_sensitive.tsv"), delimiter="\t")}
    assert got["HASLEAD"] == "sub_threshold_lead"
    assert got["ORPHAN"] == "orphan_no_pfam_similarity"
    # A curated domain settles it; the sensitive scan cannot demote it.
    assert got["CURATED"] == "has_curated_domain"


def test_sub_threshold_hits_are_never_reported_as_annotations(tmp_path):
    """The column name is load-bearing.

    Everything in a no-cut_ga scan is below Pfam's curated bar by construction.
    If it lands in a column called `pfam_domains` it will be read, and cited, as
    an annotation.
    """
    out = tmp_path / "ps"
    dom = tmp_path / "t.domtbl"
    f = ["X"] * 23
    f[0], f[3], f[12], f[13], f[22] = "HTH_38", "B|0", "0.00012", "22.5", "HTH"
    dom.write_text(" ".join(f) + "\n")
    strict = tmp_path / "s.tsv"
    _tsv(strict, ["bag", "all_domains"], [{"bag": "B", "all_domains": ""}])
    subprocess.run(
        [PY, "-I", str(ROOT / "scripts" / "pfam_sensitive_report.py"),
         "--domtbl", str(dom), "--strict", str(strict), "--out", str(out)],
        check=True, capture_output=True, text=True, env=ENV)
    hdr = open(out / "pfam_sensitive.tsv").readline().rstrip("\n").split("\t")
    assert "sensitive_leads" in hdr
    assert "pfam_domains" not in hdr
    row = next(csv.DictReader(open(out / "pfam_sensitive.tsv"), delimiter="\t"))
    assert row["curated_domains"] == ""
    assert "HTH_38" in row["sensitive_leads"]
