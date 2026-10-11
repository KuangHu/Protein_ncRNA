#!/usr/bin/env python
"""A GenBank file is read as truth, so its coordinates must be provably right.

Every other output of this pipeline is a table someone reads next to its own
documentation. A .gbk gets opened in a genome browser, lifted into a figure, and
compared against a reference annotation, with nothing to say which frame it is
in. Two things therefore have to hold and are pinned here.

The window-to-contig mapping. win_start/win_end are 1-based inclusive, and the
first version of this exporter treated win_start as 0-based. That shifted every
plus-strand feature one base left while leaving lengths, strands and the minus
strand correct -- invisible to any check except the one that round-trips the
anchor against the span written in its own id. 176 of 416 cases were wrong.

The anchor frame. Blocks are scored at x = window_pos - anchor_offset and x is
routinely negative. A feature placed at x rather than anchor_offset + x lands
thousands of bases away, or silently clamps to the window edge.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable
ENV = {"PN_ALLOW_LOGIN_NODE": "1", "PATH": "/usr/bin:/bin"}

BLOCK_COLS = ["block_id", "clade", "family", "decoy", "member_fraction",
              "n_members", "n_clade_windows", "recurrence", "coord_tightness",
              "strand_consistency", "length_consistency", "composite_score",
              "median_start", "median_len", "median_pident", "dist_to_anchor",
              "side"]
MEM_COLS = ["block_id", "clade", "decoy", "anchor_id", "start", "end", "len",
            "strand", "pident"]


def _load():
    spec = importlib.util.spec_from_file_location(
        "export_loci_gbk", ROOT / "scripts" / "export_loci_gbk.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_plus_strand_window_start_is_one_based(tmp_path):
    """The bug that put 176 of 416 plus-strand cases one base off.

    The anchor's true span is in its own id, so the mapping is checkable
    directly. A window beginning at 1-based 94601 with the anchor 5000 bases in
    puts the anchor's first base at 99601, not 99602.
    """
    g = _load().genomic
    # plus strand: anchor at window [5000, 6224) of window 94601..105824
    assert g(5000, 6224, 94601, 105824, 1) == (99601, 100824)
    # minus strand: the window is reverse-complemented, so offsets count back
    # from win_end. This branch was right the first time and must stay right.
    assert g(5000, 6224, 79543, 90766, -1) == (84543, 85766)
    # a window's full span must reproduce the window itself on both strands
    assert g(0, 11224, 94601, 105824, 1) == (94601, 105824)
    assert g(0, 11224, 79543, 90766, -1) == (79543, 90766)


def _fixture(tmp_path, *, strand=1):
    """One bag, two cases: one carrying a block upstream of the anchor, one not."""
    wd = tmp_path / "win" / "BAG1"
    wd.mkdir(parents=True)
    seq = "ACGT" * 2000                       # 8000 bp
    ao, al = 5000, 300
    ws, we = (94601, 102600) if strand == 1 else (79543, 87542)
    span = ((ws + ao, ws + ao + al - 1) if strand == 1
            else (we - (ao + al) + 1, we - ao))
    sgn = "+" if strand == 1 else "-"
    recs = []
    for k in range(2):
        recs.append({
            "anchor_id": f"GCA_00000000{k}.1|CTG{k}.1:{span[0]}-{span[1]}:{sgn}",
            "contig": f"CTG{k}.1", "win_start": ws, "win_end": we,
            "strand": strand, "anchor_offset": ao, "anchor_len": al,
            "family": "BAG1", "species": "testus_organismus",
            "window_complete": True, "seq": seq})
    (wd / "BAG1.jsonl").write_text(
        "\n".join(json.dumps(r) for r in recs) + "\n")

    def tsv(p, cols, rows):
        with open(p, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    tsv(tmp_path / "blocks.tsv", BLOCK_COLS, [{
        "block_id": "BAG1__B000", "clade": "BAG1", "family": "BAG1",
        "decoy": "False", "composite_score": "0.95", "recurrence": "0.95",
        "member_fraction": "0.95", "n_members": "1", "n_clade_windows": "2",
        "dist_to_anchor": "-400", "side": "upstream", "median_start": "-400",
        "median_len": "100"}])
    # Only the FIRST case carries the block; x is negative, as blocks upstream
    # of the anchor always are.
    tsv(tmp_path / "members.tsv", MEM_COLS, [{
        "block_id": "BAG1__B000", "clade": "BAG1", "decoy": "False",
        "anchor_id": recs[0]["anchor_id"], "start": "-400", "end": "-300",
        "len": "100", "strand": "+", "pident": "99.1"}])
    return wd, ao, al, recs


def _run(tmp_path, wd):
    out = tmp_path / "out"
    r = subprocess.run(
        [PY, "-I", str(ROOT / "scripts" / "export_loci_gbk.py"),
         "--windows", str(wd.parent), "--blocks", str(tmp_path / "blocks.tsv"),
         "--members", str(tmp_path / "members.tsv"), "--out", str(out)],
        capture_output=True, text=True, env=ENV)
    assert r.returncode == 0, r.stderr
    return out, r.stdout


def test_a_block_lands_in_the_anchor_frame_not_the_window_frame(tmp_path):
    """x = -400 is 400 bases before the anchor, not 400 before the window.

    Placing it at the raw x would put the feature off the record entirely; the
    correct position is anchor_offset + x.
    """
    from Bio import SeqIO
    wd, ao, al, recs = _fixture(tmp_path)
    out, _ = _run(tmp_path, wd)
    recs_out = {p.name: SeqIO.read(p, "genbank")
                for p in (out / "BAG1" / "gbk").glob("*.gbk")}
    assert len(recs_out) == 2
    rna = [f for r in recs_out.values() for f in r.features if f.type == "ncRNA"]
    assert len(rna) == 1, "only the case that carries the block may be annotated"
    loc = rna[0].location
    assert int(loc.start) == ao - 400 == 4600
    assert int(loc.end) == ao - 300 == 4700
    cds = [f for r in recs_out.values() for f in r.features if f.type == "CDS"]
    assert all(int(c.location.start) == ao and int(c.location.end) == ao + al
               for c in cds)


def test_a_case_without_a_block_still_gets_a_record(tmp_path):
    """Absence is a result.

    A locus carrying the protein and not the RNA is evidence about how tethered
    the association is. Emitting only the positives would turn a 99% recurrence
    and a 10% recurrence into identical-looking folders.
    """
    from Bio import SeqIO
    wd, _, _, _ = _fixture(tmp_path)
    out, _ = _run(tmp_path, wd)
    gbks = sorted((out / "BAG1" / "gbk").glob("*.gbk"))
    assert len(gbks) == 2
    empty = [p for p in gbks
             if not any(f.type == "ncRNA"
                        for f in SeqIO.read(p, "genbank").features)]
    assert len(empty) == 1
    cases = list(csv.DictReader(
        open(out / "BAG1" / "cases.tsv"), delimiter="\t"))
    assert sorted(c["n_ncrna_features"] for c in cases) == ["0", "1"]
    info = json.loads((out / "BAG1" / "bag_info.json").read_text())
    assert info["n_cases"] == 2
    assert info["n_cases_with_an_ncrna_block"] == 1
    assert info["fraction_of_cases_with_a_block"] == 0.5


def test_unverifiable_genomic_coordinates_are_omitted_not_guessed(tmp_path):
    """A wrong contig coordinate in a .gbk is worse than a missing one.

    If the mapping does not round-trip on the anchor it is not trusted for the
    blocks either, and the record says so in its COMMENT.
    """
    from Bio import SeqIO
    wd, ao, al, recs = _fixture(tmp_path)
    # Corrupt one anchor's embedded span so the round-trip must fail.
    rows = [json.loads(x) for x in
            (wd / "BAG1.jsonl").read_text().splitlines()]
    rows[0]["anchor_id"] = "GCA_000000000.1|CTG0.1:1-2:+"
    (wd / "BAG1.jsonl").write_text(
        "\n".join(json.dumps(r) for r in rows) + "\n")
    # members key on anchor_id, so re-point the member at the corrupted one
    txt = (tmp_path / "members.tsv").read_text().replace(
        recs[0]["anchor_id"], "GCA_000000000.1|CTG0.1:1-2:+")
    (tmp_path / "members.tsv").write_text(txt)

    out, stdout = _run(tmp_path, wd)
    assert "1 without verified genomic coordinates" in stdout
    bad = [p for p in (out / "BAG1" / "gbk").glob("*.gbk")
           if "OMITTED" in SeqIO.read(p, "genbank").annotations["comment"]]
    assert len(bad) == 1
    r = SeqIO.read(bad[0], "genbank")
    # "contig strand +" may appear; a coordinate span "contig 123..456" may not.
    import re as _re
    for feat in r.features:
        if feat.type in ("CDS", "ncRNA"):
            assert not _re.search(r"contig \d+\.\.\d+",
                                  feat.qualifiers["note"][0])
    cases = {c["anchor_id"]: c for c in csv.DictReader(
        open(out / "BAG1" / "cases.tsv"), delimiter="\t")}
    assert cases["GCA_000000000.1|CTG0.1:1-2:+"][
        "genomic_coords_verified"] == "False"


def test_minus_strand_cases_verify_too(tmp_path):
    """The minus branch counts back from win_end and must not regress."""
    wd, _, _, _ = _fixture(tmp_path, strand=-1)
    out, stdout = _run(tmp_path, wd)
    assert "0 without verified genomic coordinates" in stdout


def test_shard_merge_ignores_io_paths_but_not_scoring_params(tmp_path):
    """Sharding changes every shard's paths and none of its scoring.

    Comparing whole params dicts fires on every shard of every correct run,
    which trains the reader to ignore the one check that would catch a shard
    scored at the wrong threshold -- a merged table that silently mixes two
    min_recurrence settings and whose floor is then computed across both.
    """
    import json as _json
    m = importlib.util.spec_from_file_location(
        "merge_block_shards", ROOT / "scripts" / "merge_block_shards.py")
    mod = importlib.util.module_from_spec(m)
    m.loader.exec_module(mod)

    base = {"windows": "/a/shard_00/real", "decoys": "/a/shard_00/decoys",
            "out": "/a/shard_00/blocks", "min_recurrence": "0.5",
            "references": "3"}
    other = dict(base, windows="/a/shard_01/real",
                 decoys="/a/shard_01/decoys", out="/a/shard_01/blocks")
    assert mod.scoring_params(base) == mod.scoring_params(other)

    wrong = dict(other, min_recurrence="0.25")
    assert mod.scoring_params(base) != mod.scoring_params(wrong)

    # end to end: two shards differing only in paths must merge cleanly
    for i, prm in enumerate((base, other)):
        d = tmp_path / f"shard_{i:02d}"
        (d / "blocks").mkdir(parents=True)
        for t in ("blocks.tsv", "block_members.tsv", "clade_stats.tsv",
                  "reference_recovery.tsv"):
            (d / t).write_text("col\n" + f"v{i}\n")
        (d / "discovery_summary.json").write_text(_json.dumps({
            "real_clades": 1, "decoy_clades": 1, "real_blocks": 0,
            "decoy_blocks": 0, "complete": True, "failures": [],
            "params": prm}))
    r = subprocess.run(
        [PY, "-I", str(ROOT / "scripts" / "merge_block_shards.py"),
         "--shards", str(tmp_path / "shard_00"), str(tmp_path / "shard_01"),
         "--out", str(tmp_path / "merged"), "--expect-clades", "4"],
        capture_output=True, text=True, env=ENV)
    assert r.returncode == 0, r.stderr + r.stdout
    got = _json.loads((tmp_path / "merged" / "discovery_summary.json").read_text())
    assert got["complete"] is True
    assert got["merge_problems"] == []
