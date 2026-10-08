#!/usr/bin/env python
"""Block-caller coordinate arithmetic.

The step that is most likely to be silently wrong here is mapping a called
block from the reference axis back onto each member, because a minus-strand
BLAST HSP relates the two axes in opposite directions while both intervals are
stored ascending. A sign error there returns a slice of the right *length* at
the wrong *place*, which no summary statistic would show: recurrence, block
length and member count all stay the same, and only the emitted sequence is
wrong. So the mirror property is pinned by a test rather than by reading BLAST
output.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

_spec = importlib.util.spec_from_file_location(
    "_db", ROOT / "scripts" / "discover_blocks.py")
_db = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_db)

trim_member, subtract = _db.trim_member, _db.subtract


def test_plus_strand_trim_is_left_to_left():
    # Reference segment [100, 200) is member segment [500, 600). The block
    # covers [120, 180), so 20 is cut from the left and 20 from the right.
    assert trim_member(100, 200, 120, 180, 500, 600, "+") == (520, 580)


def test_minus_strand_trim_is_mirrored():
    # Same HSP, opposite orientation: reference-left answers to member-right.
    assert trim_member(100, 200, 120, 180, 500, 600, "-") == (520, 580)
    # An asymmetric block is where the two differ. Cutting 30 from the
    # reference's left end must cut 30 from the *member's right* end.
    assert trim_member(100, 200, 130, 200, 500, 600, "+") == (530, 600)
    assert trim_member(100, 200, 130, 200, 500, 600, "-") == (500, 570)


def test_strands_are_exact_mirrors_of_each_other():
    """For every sub-block, the plus and minus results reflect through the HSP.

    Formally: minus_start - ma == mb - plus_end, and minus_end - ma ==
    mb - plus_start. This is the property that fails if either branch is
    written with the wrong sign.
    """
    a, b, ma, mb = 100, 400, 7000, 7300
    for ba in range(100, 400, 17):
        for bb in range(ba + 10, 401, 23):
            pa, pb = trim_member(a, b, ba, bb, ma, mb, "+")
            na, nb = trim_member(a, b, ba, bb, ma, mb, "-")
            assert nb - ma == mb - pa, (ba, bb, pa, pb, na, nb)
            assert na - ma == mb - pb, (ba, bb, pa, pb, na, nb)
            # Length is orientation-independent; only the position moves.
            assert pb - pa == nb - na
            # Both stay inside the HSP.
            assert ma <= pa < pb <= mb and ma <= na < nb <= mb


def test_block_wider_than_segment_is_not_extended():
    """A block overhanging the segment trims to the segment, never past it."""
    for strand in ("+", "-"):
        assert trim_member(100, 200, 50, 300, 500, 600, strand) == (500, 600)


def test_subtract_trims_rather_than_discards():
    """An HSP spanning a gene and the intergenic stretch beside it keeps the
    intergenic part. Dropping such an HSP whole is what reduced a 493-member
    clade to 154 usable members when these windows are ~86% coding."""
    assert subtract(0, 100, [(40, 60)]) == [(0, 40), (60, 100)]
    assert subtract(0, 100, [(0, 100)]) == []
    assert subtract(0, 100, []) == [(0, 100)]


def test_partial_step4_run_is_refused_downstream():
    """A step-4 run that lost clades must not become an FDR denominator.

    The failure this guards is quiet by construction: `real_clades` and
    `decoy_clades` in discovery_summary.json count *attempted* clades, so a run
    that crashed on 20 decoys still reports 129 and the empirical FDR simply
    improves. The contract is therefore a `complete` flag written by
    discover_blocks.py and checked by candidate_table.py.
    """
    src = (ROOT / "scripts" / "discover_blocks.py").read_text()
    assert '"complete": not fails' in src, \
        "discover_blocks.py must record whether every clade succeeded"
    assert "failure_report.tsv" in src and "allow_partial" in src, \
        "a failed clade must produce a report and a non-zero exit by default"
    # The refusal has to be the default, not an opt-in.
    i = src.index("if not args.allow_partial")
    assert "raise SystemExit" in src[i:i + 400], \
        "partial step-4 output must exit non-zero unless --allow-partial"

    ct = (ROOT / "scripts" / "candidate_table.py").read_text()
    assert 'd.get("complete") is False' in ct and "raise SystemExit" in ct, \
        "candidate_table.py must refuse an incomplete step-4 run"


if __name__ == "__main__":
    import traceback
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"  ok   {name}")
            except Exception:
                fails += 1
                print(f"  FAIL {name}")
                traceback.print_exc()
    print("FAILED" if fails else "all block coordinate tests passed")
    raise SystemExit(1 if fails else 0)
