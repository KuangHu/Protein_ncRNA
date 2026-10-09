#!/usr/bin/env python
"""The spike-in lane must not leak into the discovery result.

IS110, group II introns, IS200/IS605 and TnpB/Cas12f are insertion systems that
are already known to carry a structured RNA. In this run they are controls: they
measure whether the blind caller recovers a signal known to be present. If one
of them reaches the discovery table it is reported as a novel protein-ncRNA
association when it is nothing of the kind, and the headline of the whole run is
wrong.

Two separate places can leak:
  * stratification -- a bag carrying BOTH a known-RNA marker and a generic
    transposase domain must land in the control lane, not the discovery lane.
  * reporting -- the discovery table must filter on lane, not on score.

Also pinned here: the score floor is recalibrated from this run's own decoys,
and genomic_background is excluded from that calibration. It produced zero
blocks in all 18 pilot cells, so including it would pull any percentile toward
zero and manufacture a floor nothing has to clear.
"""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

PY = sys.executable
ENV = {"PN_ALLOW_LOGIN_NODE": "1", "PATH": "/usr/bin:/bin"}


def _load(name: str):
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_known_rna_marker_beats_a_generic_transposase_domain():
    """The precedence that keeps a control out of the discovery table.

    A group II intron RT sitting next to a DDE domain, or an IS110 bag that also
    hits `rve`, must be a control. Sorting on whichever domain happened to score
    best would put it in whichever lane the HMM gods chose.
    """
    s = _load("stratify_bags")
    assert s.assign(["RVT_1", "DDE_Tnp_1", "rve"])[0] == "known_rna_mge"
    assert s.assign(["DDE_Tnp_1", "DEDD_Tnp_IS110"])[0] == "known_rna_mge"
    assert s.assign(["Y1_Tnp", "HTH_21"])[0] == "known_rna_mge"
    assert s.assign(["Cas12f1-like_TNB", "rve_3"])[0] == "known_rna_mge"
    # ... and only then does a generic transposase get the discovery lane.
    assert s.assign(["DDE_Tnp_1", "HTH_21"]) == ("discovery", "DDE_transposase")


def test_a_bag_with_only_a_generic_domain_is_unknown_not_a_transposase():
    """An HTH says the protein binds DNA, not what it is.

    Calling such a bag a transposase would put it in the discovery lane with a
    label it did not earn; calling it host-like would hide the most novel thing
    the run can return. It is `unknown`, and it is its own lane.
    """
    s = _load("stratify_bags")
    assert s.assign(["HTH_21", "HTH_38"]) == ("unknown", "no_specific_domain")
    assert s.assign([]) == ("unknown", "no_specific_domain")
    assert s.assign(["MCPsignal"])[0] == "host_like"


def test_multiple_known_markers_are_all_named():
    s = _load("stratify_bags")
    lane, grp = s.assign(["DEDD_Tnp_IS110", "GIIM"])
    assert lane == "known_rna_mge"
    assert grp == "IS110+group_II_intron"


def _blocks(tmp: Path, rows: list[tuple[str, str, float, float]]) -> Path:
    p = tmp / "blocks.tsv"
    with open(p, "w") as fh:
        fh.write("block_id\tclade\tfamily\tdecoy\tcomposite_score\t"
                 "member_fraction\n")
        for i, (clade, decoy, score, mf) in enumerate(rows):
            fh.write(f"B{i:03d}\t{clade}\t{clade}\t{decoy}\t{score}\t{mf}\n")
    return p


def test_floor_ignores_genomic_background(tmp_path):
    """A null that never produces a block cannot set a floor.

    genomic_background decoys scored 0 blocks in all 18 pilot cells. Were they
    admitted to the calibration, a lone low-scoring one would drag the
    percentile down and every real block would clear a floor nothing defended.
    """
    blocks = _blocks(tmp_path, [
        ("CDS0001", "False", 0.95, 0.9),
        ("DECOY_genomic_background__r0__CDS0001", "True", 0.01, 0.6),
        ("DECOY_within_bag_shuffle__r0__CDS0001", "True", 0.42, 0.6),
        ("DECOY_same_family_permutation__r0__CDS0001", "True", 0.61, 0.7),
    ])
    out = tmp_path / "floor"
    r = subprocess.run(
        [PY, str(ROOT / "scripts" / "arm2_floor.py"),
         "--blocks", str(blocks), "--out", str(out)],
        capture_output=True, text=True, env=ENV)
    assert r.returncode == 0, r.stderr
    fl = json.loads((out / "score_floor.json").read_text())
    # The strict nulls reached 0.61; the 0.01 background block must not appear
    # anywhere in the calibration.
    assert fl["floor_no_decoy_support"] == 0.61
    assert fl["n_strict_decoy_blocks"] == 2
    assert fl["per_mode"]["genomic_background"]["n_blocks"] == 1
    assert fl["n_real_blocks_above_floor_no_decoy_support"] == 1
    # The legacy floor is recorded but never used as the answer.
    assert fl["legacy_genomic_floor"] == 0.80


def test_no_strict_decoy_block_means_no_floor_and_a_warning(tmp_path):
    blocks = _blocks(tmp_path, [
        ("CDS0001", "False", 0.95, 0.9),
        ("DECOY_genomic_background__r0__CDS0001", "True", 0.01, 0.6),
    ])
    out = tmp_path / "floor"
    subprocess.run([PY, str(ROOT / "scripts" / "arm2_floor.py"),
                    "--blocks", str(blocks), "--out", str(out)],
                   capture_output=True, text=True, env=ENV, check=True)
    fl = json.loads((out / "score_floor.json").read_text())
    assert fl["floor_no_decoy_support"] is None
    assert "warning" in fl


def _lane_fixture(tmp: Path):
    """One control bag and one discovery bag, control scoring higher."""
    (tmp / "decoys").mkdir()
    for bag in ("CDS0001", "CDS0002"):
        for m in ("within_bag_shuffle", "same_family_permutation"):
            (tmp / "decoys" / f"DECOY_{m}__rep0__{bag}.fasta").write_text(">x\nA\n")
    blocks = tmp / "blocks"
    blocks.mkdir()
    _blocks(blocks, [
        ("CDS0001", "False", 0.99, 0.95),            # known_rna_mge control
        ("CDS0002", "False", 0.88, 0.80),            # discovery
        ("DECOY_within_bag_shuffle__rep0__CDS0001", "True", 0.10, 0.3),
        ("DECOY_within_bag_shuffle__rep0__CDS0002", "True", 0.12, 0.3),
        ("DECOY_same_family_permutation__rep0__CDS0001", "True", 0.20, 0.4),
        ("DECOY_same_family_permutation__rep0__CDS0002", "True", 0.25, 0.4),
    ])
    with open(blocks / "clade_stats.tsv", "w") as fh:
        fh.write("clade\tdecoy\n")
        for b in ("CDS0001", "CDS0002"):
            fh.write(f"{b}\tFalse\n")
    with open(tmp / "lanes.tsv", "w") as fh:
        fh.write("bag\tlane\tlane_group\tn_distinct\tfamily_size\n")
        fh.write("CDS0001\tknown_rna_mge\tIS110\t500\t6\n")
        fh.write("CDS0002\tdiscovery\tDDE_transposase\t400\t6\n")
    with open(tmp / "strata.tsv", "w") as fh:
        fh.write("bag\tlane\tlane_group\ttop_domain\ttop_domain_frac_of_reps\t"
                 "all_domains\n")
        fh.write("CDS0001\tknown_rna_mge\tIS110\tDEDD_Tnp_IS110\t1.0\t"
                 "DEDD_Tnp_IS110(3)\n")
        fh.write("CDS0002\tdiscovery\tDDE_transposase\tDDE_Tnp_1\t1.0\t"
                 "DDE_Tnp_1(3)\n")
    (tmp / "floor.json").write_text(json.dumps({
        "floor_no_decoy_support": 0.25,
        "floor_backed_by_n_decoy_blocks": 4,
        "legacy_genomic_floor": 0.80}))
    return blocks


def test_controls_never_reach_the_discovery_table(tmp_path):
    """The highest-scoring bag in the run is a spike-in. It must not be the
    headline, and it must still be reported -- in its own table."""
    blocks = _lane_fixture(tmp_path)
    out = tmp_path / "report"
    r = subprocess.run(
        [PY, str(ROOT / "scripts" / "arm2_lane_report.py"),
         "--blocks", str(blocks), "--decoys", str(tmp_path / "decoys"),
         "--lanes", str(tmp_path / "lanes.tsv"),
         "--strata", str(tmp_path / "strata.tsv"),
         "--floor", str(tmp_path / "floor.json"), "--out", str(out)],
        capture_output=True, text=True, env=ENV)
    assert r.returncode == 0, r.stderr
    disc = (out / "discovery_candidates.tsv").read_text()
    assert "CDS0002" in disc
    assert "CDS0001" not in disc, (
        "the IS110 control reached the discovery table; it scored highest, "
        "which is exactly when the filter matters")
    # It is not dropped, only relocated.
    assert "CDS0001" in (out / "positive_control_recovery.tsv").read_text()
    assert "CDS0001" in (out / "all_bags_annotated.tsv").read_text()
    # Recovery is reported for the spike-in lane.
    assert "positive control recovery" in r.stdout


def test_decoy_status_separates_an_absent_null_from_a_clean_one(tmp_path):
    """`not_available` and `clean` are opposite findings and must not share a
    cell. Arm 1 labelled 1,497 blockless decoy clades `strong`; that is the bug
    this column exists to make impossible."""
    blocks = _lane_fixture(tmp_path)
    # Remove CDS0002's permutation decoy so the null is genuinely unavailable.
    (tmp_path / "decoys" /
     "DECOY_same_family_permutation__rep0__CDS0002.fasta").unlink()
    out = tmp_path / "report2"
    subprocess.run(
        [PY, str(ROOT / "scripts" / "arm2_lane_report.py"),
         "--blocks", str(blocks), "--decoys", str(tmp_path / "decoys"),
         "--lanes", str(tmp_path / "lanes.tsv"),
         "--strata", str(tmp_path / "strata.tsv"),
         "--floor", str(tmp_path / "floor.json"), "--out", str(out)],
        capture_output=True, text=True, env=ENV, check=True)
    rows = [l.split("\t") for l in
            (out / "all_bags_annotated.tsv").read_text().splitlines()]
    hdr = rows[0]
    col = hdr.index("same_family_permutation__status")
    got = {r[0]: r[col] for r in rows[1:]}
    assert got["CDS0002"] == "not_available"
    assert got["CDS0001"] == "clean"


def _win(root: Path, bag: str, n: int = 2) -> None:
    d = root / bag
    (d / "proteins").mkdir(parents=True)
    (d / f"{bag}.fasta").write_text(
        "".join(f">w{i} family={bag}\nACGT\n" for i in range(n)))
    (d / f"{bag}.jsonl").write_text(
        "".join(json.dumps({"anchor_id": f"w{i}", "family": bag}) + "\n"
                for i in range(n)))
    (d / "proteins" / f"{bag}.faa").write_text(f">w0\nMKV\n")


def test_a_pool_bag_is_staged_under_every_family_it_belongs_to(tmp_path):
    """Remote-homology pools overlap, and the overlap must survive staging.

    46 of 136 pool bags are homologs of more than one target; CDS02588 serves
    five. A bag -> family map assigns each such bag to one family and silently
    deletes it from the others, shrinking pools that were sized and labelled on
    the full membership -- so a pool reported `strong` would be scored as
    something smaller, with no error anywhere.
    """
    wins = tmp_path / "w"
    wins.mkdir()
    for b in ("TGT1", "TGT2", "SHARED", "ONLY1", "ONLY2"):
        _win(wins, b)
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text(
        "family\tbag_id\trole\tfamily_size\n"
        "H000\tTGT1\ttarget\t3\nH000\tSHARED\tpool\t3\nH000\tONLY1\tpool\t3\n"
        "H001\tTGT2\ttarget\t3\nH001\tSHARED\tpool\t3\nH001\tONLY2\tpool\t3\n")
    subj = tmp_path / "subj.txt"
    subj.write_text("TGT1\nTGT2\n")
    out = tmp_path / "staged"
    r = subprocess.run(
        [PY, str(ROOT / "scripts" / "stage_family_windows.py"),
         "--windows", str(wins), "--pairs", str(pairs),
         "--subjects", str(subj), "--out", str(out), "--min-size", "2"],
        capture_output=True, text=True, env=ENV)
    assert r.returncode == 0, r.stderr
    staged = {p.stem for p in out.glob("*.fasta")}
    assert "H000__SHARED" in staged and "H001__SHARED" in staged, (
        f"the shared pool bag reached only one family: {sorted(staged)}")
    # Both pools keep their full size, which is what they were labelled on.
    for fam in ("H000", "H001"):
        assert sum(1 for s in staged if s.startswith(fam + "__")) == 3
    # The family= tag is rewritten per staged copy, not left as the bag name.
    assert "family=H000" in (out / "H000__SHARED.fasta").read_text()
    assert "family=H001" in (out / "H001__SHARED.fasta").read_text()
