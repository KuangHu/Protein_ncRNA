#!/usr/bin/env python
"""The permutation null must actually break relatedness.

A decoy clade is only a null if its members are *not* each other's relatives.
Drawing uniformly from the family pool does not guarantee that: the pool is
dominated by its largest clades, so a decoy can be mostly one real clade under a
decoy name. Measured on the 43-clade set before this was fixed, Cas3 decoys drew
67-71% of their members from a single source clade, and Cas3 was also the one
family whose "decoys" cleared the score floor -- which is what a null made of
real clades would do.

The bias runs in the conservative direction, which is why it stayed hidden: a
decoy that keeps relatedness produces *more* blocks and so *raises* the measured
FDR. A family reporting zero decoy blocks is therefore safe; a family reporting
many may just have had a degenerate null. Hence the diagnostics, and hence
labelling rather than silently reporting.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SRC = (ROOT / "scripts" / "permutation_null.py").read_text()


def test_decoy_excludes_the_clade_it_mirrors_by_default():
    assert "allow_self_source or s != stem" in SRC, \
        "a decoy must not draw from the clade it mirrors unless asked to"
    assert '"--allow-self-source"' in SRC


def test_per_source_quota_exists():
    """Without a quota, round-robin or not, one big clade can fill the decoy."""
    assert "quota = -(-n // len(srcs))" in SRC, "missing ceil(N / n_sources) quota"
    assert "taken[s] >= quota" in SRC, "the quota must actually gate the draw"


def test_degeneracy_is_recorded_not_hidden():
    for k in ("n_source_clades", "max_source_fraction", "null_strength",
              "degenerate"):
        assert k in SRC, f"decoy manifest must record {k}"
    assert "min_source_clades" in SRC and "max_source_fraction" in SRC


def test_candidate_table_labels_null_strength():
    ct = (ROOT / "scripts" / "candidate_table.py").read_text()
    assert "null_strength" in ct, \
        "the per-family FDR table must carry its null's strength"
    assert "families_with_degenerate_null" in ct


def test_rscape_dir_is_cleared_before_each_run():
    """The power file is located by glob, so a stale one would be read as this
    run's power analysis -- and power is what separates `underpowered` from
    `tested_no_covariation`."""
    fb = (ROOT / "scripts" / "fold_blocks.py").read_text()
    i = fb.index('d = outdir / "rscape"')
    seg = fb[i:i + 500]
    assert "shutil.rmtree(d)" in seg, "rscape outdir must be wiped before reuse"
    assert "for pf in sorted(d.glob" not in fb, \
        "do not take an arbitrary .power file; wipe the dir and warn on >1"
    assert "len(powers) > 1" in fb


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
    print("FAILED" if fails else "all null tests passed")
    raise SystemExit(1 if fails else 0)


def _arm2():
    spec = importlib.util.spec_from_file_location(
        "_a2", ROOT / "scripts" / "arm2_nulls.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_shuffled_pseudo_anchor_clears_the_real_one_in_one_frame():
    """Both offsets must be reported in the frame the sequence ends up in.

    shuffle_anchor randomises orientation. Flipping the pseudo-anchor's offset
    while leaving `real_anchor_offset` in the pre-flip frame leaves half the
    records describing two different coordinate systems, so `|pseudo - real|`
    silently stops meaning separation -- which is the one property the >=2 kb
    rule exists to guarantee.
    """
    import random

    a2 = _arm2()
    alen, wlen, sep = 600, 11000, 2000
    recs = [{"anchor_id": f"a{i}", "seq": "ACGT" * (wlen // 4),
             "anchor_offset": 5000, "anchor_len": alen, "family": "F",
             "species": ""} for i in range(40)]
    out = a2.shuffle_anchor(recs, random.Random(0), 0, sep)
    assert out, "no shuffled decoys produced"
    flipped = 0
    for r in out:
        p, q = r["anchor_offset"], r["real_anchor_offset"]
        assert 0 <= p <= len(r["seq"]) - alen
        assert 0 <= q <= len(r["seq"]) - alen
        # Separation holds between the two ends, not just the two starts.
        assert p + alen + sep <= q or q + alen + sep <= p, (p, q)
        if q != r["real_anchor_offset_input_frame"]:
            flipped += 1
    assert flipped, "orientation was never randomised in 40 draws"


def test_shuffled_decoys_carry_their_anchor_proteins(tmp_path):
    """Every shuffled anchor id must resolve to a protein under its new name.

    shuffle_anchor renames each anchor `<id>__sh<rep>`. discover_blocks looks up
    anchor proteins by exact id and silently falls back to window DNA when the
    lookup misses, so without re-emitting the proteins under the shuffled ids
    this arm would pick references in a different space from the real arm --
    and it is the arm that is supposed to differ from real in exactly one way,
    the anchor frame.
    """
    import random

    a2 = _arm2()
    alen, wlen = 600, 11000
    recs = [{"anchor_id": f"a{i}", "seq": "ACGT" * (wlen // 4),
             "anchor_offset": 5000, "anchor_len": alen, "family": "F",
             "species": ""} for i in range(12)]
    pin, pout = tmp_path / "in", tmp_path / "out"
    pin.mkdir()
    with open(pin / "BAG.faa", "w") as fh:
        for i in range(12):
            fh.write(f">a{i}\nM{'ACDEFGHIKL'[i % 10] * 60}\n")

    d = a2.shuffle_anchor(recs, random.Random(3), 0, 2000)
    n = a2.write_shuffled_proteins(pin, pout, "BAG", d, 0)
    assert n == len(d), f"{n} proteins for {len(d)} shuffled windows"

    emitted = {line[1:].split()[0]
               for line in open(pout / "DECOY_within_bag_shuffle__rep0__BAG.faa")
               if line.startswith(">")}
    assert emitted == {r["anchor_id"] for r in d}
    assert all(i.endswith("__sh0") for i in emitted)
