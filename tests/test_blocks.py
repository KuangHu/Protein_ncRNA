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


def test_gapped_hsp_trim_never_crosses_coordinates():
    """A member trimmed to nothing is dropped, not returned inverted.

    The trim is measured on the reference axis and applied to the member axis,
    but a BLAST HSP is gapped, so those two axes are different lengths. Where
    the member carries a deletion and the overlap with the block is small, the
    left and right trims together exceed the member's own span and the result
    comes back with `end < start` -- an interval denoting no sequence at all.

    It reached production. 85 of 66,603 members in the published genomic run
    were inverted (73 of them on the PLUS strand, so this is not the
    minus-strand mirror bug), and on the insertion corpus one block came back
    with median_len = -67, which drove `length_consistency = 1 - mad/med_len`
    above 1 and produced a composite score of 2.15 on a scale ending at 1.
    """
    # Scaling is the fix; dropping is the backstop. A 100 bp reference segment
    # against a 20 bp member is a heavily gapped HSP, and trimming 60 bp off the
    # reference axis verbatim would cut more than the member has -- scaled, it
    # cuts 12 and the interval stays valid.
    for strand in ("+", "-"):
        t = trim_member(100, 200, 160, 170, 500, 520, strand)
        assert t is not None and t[1] > t[0], (strand, t)

    # Where even the scaled trim leaves nothing, the member is dropped rather
    # than returned crossed.
    for strand in ("+", "-"):
        assert trim_member(100, 200, 160, 170, 500, 501, strand) is None, strand

    # And over a sweep, nothing may ever come back crossed.
    bad = []
    for ref_span, mem_span in ((100, 10), (100, 40), (100, 250), (50, 50)):
        a, b, ma, mb = 100, 100 + ref_span, 500, 500 + mem_span
        for ba in range(a - 20, b + 20, 7):
            for bb in range(ba + 1, b + 40, 11):
                if min(b, bb) - max(a, ba) <= 0:
                    continue
                for strand in ("+", "-"):
                    t = trim_member(a, b, ba, bb, ma, mb, strand)
                    if t is None:
                        continue
                    fa, fb = t
                    if fb <= fa:
                        bad.append((ref_span, mem_span, ba, bb, strand, fa, fb))
                    if fa < ma or fb > mb:
                        bad.append(("escaped HSP", ref_span, mem_span, ba, bb,
                                    strand, fa, fb))
    assert not bad, f"{len(bad)} crossed/escaped intervals, e.g. {bad[:3]}"


def test_trim_scales_with_the_hsp_length_ratio():
    """The member axis is cut in proportion, not by the reference's bp count."""
    # Member span is half the reference span, so a 40 bp reference trim is a
    # 20 bp member trim.
    assert trim_member(0, 100, 40, 100, 0, 50, "+") == (20, 50)
    # Member span is double, so the same 40 bp reference trim cuts 80 bp.
    assert trim_member(0, 100, 40, 100, 0, 200, "+") == (80, 200)
    # Ungapped (1:1) keeps the original behaviour exactly.
    assert trim_member(100, 200, 120, 180, 500, 600, "+") == (520, 580)


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


def test_identity_estimate_deduplicates_before_sampling():
    """Median pairwise identity must measure divergence, not copy number.

    These clades are mostly duplicates -- one GroupII_RT clade has 799 members
    and 29 distinct proteins. Sampling raw members therefore draws almost only
    duplicate-vs-duplicate pairs and the median goes to 1.0 whatever the clade
    actually looks like. The same clade measured 0.451 from an ordered sample
    and 1.000 from a random one; the second crossed the `identity > 0.9`
    rejection gate and cost GroupII_RT 3,415 alignable windows.
    """
    src = (ROOT / "scripts" / "clade_decompose.py").read_text()
    i = src.index("def approx_median_identity")
    body = src[i:i + 3000]
    assert "first.setdefault(seqs[m], m)" in body, \
        "the identity sample must be drawn from distinct sequences"
    assert "hash_sample(uniq" in body, \
        "sample the deduplicated set, not the raw member list"
    assert "members[:cap]" not in body, "no order-dependent slice"


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
