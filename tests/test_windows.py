#!/usr/bin/env python
"""Strand normalization is the step most likely to be silently wrong, so it gets
an explicit test: the same locus encoded on either strand must produce the same
oriented window, with a marker block landing at the same offset in both.
"""

from __future__ import annotations

import csv
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_ncrna.anchors import Anchor
from protein_ncrna.orfs import ORF
from protein_ncrna.seqio import revcomp
from protein_ncrna.windows import extract_window

CFG = {"window": {"flank_bp": 100, "orient_to_anchor_strand": True, "min_contig_bp": 0}}
MARKER = "ACGTACGTAAATTTCCCGGG"


def _anchor(start, end, strand, contig="c1"):
    return Anchor(orf=ORF(contig, start, end, strand, 0, "M" * 120),
                  family="RT", model="RVT_1", evalue=1e-20, score=200.0)


def test_strand_symmetry():
    rng = random.Random(7)
    bg = "".join(rng.choice("ACGT") for _ in range(600))

    # Marker 60 bp upstream of a 150 bp anchor that starts at 201 on the + strand.
    a_start, a_len = 201, 150
    m_start = a_start - 60 - len(MARKER)
    fwd = bg[:m_start - 1] + MARKER + bg[m_start - 1 + len(MARKER):]
    a_end = a_start + a_len - 1

    w_plus = extract_window(_anchor(a_start, a_end, +1), fwd, "g1", CFG)

    # Same locus on the opposite strand: revcomp the contig and remap coordinates.
    rc = revcomp(fwd)
    n = len(fwd)
    r_start, r_end = n - a_end + 1, n - a_start + 1
    w_minus = extract_window(_anchor(r_start, r_end, -1), rc, "g2", CFG)

    assert w_plus.seq == w_minus.seq, "oriented windows differ between strands"
    assert w_plus.anchor_offset == w_minus.anchor_offset, "anchor offset differs"
    assert w_plus.left_flank_bp == w_minus.left_flank_bp
    assert w_plus.right_flank_bp == w_minus.right_flank_bp

    # The anchor must actually sit where anchor_offset says it does.
    assert w_plus.anchor_len == a_len
    # The marker must sit at the same oriented offset in both.
    assert w_plus.seq.find(MARKER) == w_minus.seq.find(MARKER) >= 0
    gap = w_plus.anchor_offset - (w_plus.seq.find(MARKER) + len(MARKER))
    assert gap == 60, f"marker-to-anchor gap {gap}, expected 60"
    print("ok  strand symmetry: anchor_offset="
          f"{w_plus.anchor_offset}, marker gap={gap}")


def test_contig_edge_truncation():
    rng = random.Random(11)
    bg = "".join(rng.choice("ACGT") for _ in range(300))
    # Anchor 20 bp from the contig start, so the left flank cannot be filled.
    w = extract_window(_anchor(21, 170, +1), bg, "g1", CFG)
    assert w.win_start == 1 and w.left_flank_bp == 20
    assert w.truncated_left and not w.truncated_right
    assert w.anchor_offset == 20
    assert len(w.seq) == w.win_end - w.win_start + 1

    # On the minus strand the truncation must flip to the other side.
    rc = revcomp(bg)
    n = len(bg)
    w2 = extract_window(_anchor(n - 170 + 1, n - 21 + 1, -1), rc, "g2", CFG)
    assert w2.seq == w.seq and w2.truncated_left and not w2.truncated_right
    print("ok  contig-edge truncation flips with strand")


def test_min_contig_filter():
    cfg = {"window": dict(CFG["window"], min_contig_bp=1000)}
    assert extract_window(_anchor(21, 170, +1), "A" * 300, "g1", cfg) is None
    print("ok  short contigs dropped")


if __name__ == "__main__":
    test_strand_symmetry()
    test_contig_edge_truncation()
    test_min_contig_filter()
    print("all window tests passed")


def test_cap_diagnostics_reads_both_protein_layouts(tmp_path):
    """nr100 must be computable on either window producer's protein layout.

    extract_windows.py writes `proteins/<family>.faa`; bags_to_windows.py writes
    `proteins/<clade>.faa`. cap_diagnostics.py keyed only on the family name, so
    on the bag corpus it loaded nothing and reported every clade as retaining 0
    of 0 distinct proteins -- a table that looks exactly like a cap that cost
    nothing, which is the one conclusion it exists to rule out.
    """
    import os
    import subprocess
    import sys as _sys

    root = Path(__file__).resolve().parents[1]
    # Both fixtures are a handful of 41-aa sequences; the SLURM guard exists to
    # stop real corpora running on a login node, and this is what its documented
    # override is for.
    env = {**os.environ, "PN_ALLOW_LOGIN_NODE": "1"}
    for layout, name in (("family", "FAM"), ("clade", "FAM__C1")):
        w = tmp_path / layout
        (w / "proteins").mkdir(parents=True)
        with open(w / "cap_manifest.tsv", "w") as fh:
            fh.write("clade_id\tanchor_id\tselected\treason\n")
            for i in range(6):
                fh.write(f"FAM__C1\ta{i}\t{str(i < 3)}\tstratified_kmer\n")
        with open(w / "proteins" / f"{name}.faa", "w") as fh:
            for i in range(6):
                fh.write(f">a{i}\n{'MA' * 20}{'CDEFGH'[i]}\n")
        r = subprocess.run(
            [_sys.executable, str(root / "scripts" / "cap_diagnostics.py"),
             "--windows", str(w), "--census", str(w), "--threads", "1",
             "--out", str(w / "cap_diagnostics.tsv")],
            capture_output=True, text=True, env=env)
        assert r.returncode == 0, f"{layout} layout failed: {r.stderr[-2000:]}"
        rows = list(csv.DictReader(open(w / "cap_diagnostics.tsv"),
                                   delimiter="\t"))
        assert len(rows) == 1
        assert int(rows[0]["pre_cap_nr100"]) == 6, rows[0]
        assert int(rows[0]["post_cap_nr100"]) == 3, rows[0]


def test_cap_diagnostics_refuses_to_report_zero_proteins(tmp_path):
    """With no protein file at all it must fail, not print nr100 0 -> 0."""
    import os
    import subprocess
    import sys as _sys

    root = Path(__file__).resolve().parents[1]
    # Both fixtures are a handful of 41-aa sequences; the SLURM guard exists to
    # stop real corpora running on a login node, and this is what its documented
    # override is for.
    env = {**os.environ, "PN_ALLOW_LOGIN_NODE": "1"}
    w = tmp_path / "empty"
    (w / "proteins").mkdir(parents=True)
    with open(w / "cap_manifest.tsv", "w") as fh:
        fh.write("clade_id\tanchor_id\tselected\treason\n")
        fh.write("FAM__C1\ta0\tTrue\tunder_cap\n")
    r = subprocess.run(
        [_sys.executable, str(root / "scripts" / "cap_diagnostics.py"),
         "--windows", str(w), "--census", str(w), "--threads", "1",
         "--out", str(w / "cap_diagnostics.tsv")],
        capture_output=True, text=True, env=env)
    assert r.returncode != 0
    assert "no anchor proteins loaded" in r.stderr
