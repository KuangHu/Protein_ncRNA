#!/usr/bin/env python
"""Subset selection must not depend on input order.

Every cap in this pipeline used to be a list slice: the first 500 anchors per
clade, the first 300 members for an identity estimate, the first 40 for a
tie-break. A slice is not a random sample -- it follows manifest, accession and
submission order, which correlate with the subclade structure the cap is
supposed to be spreading across. The failure is silent twice over: the windows
that are never extracted leave no trace, and when R-scape later calls a block
underpowered there is no way to separate a biologically tight clade from a
narrowly sampled one.

Two different contracts are tested, because the two samplers exist for opposite
reasons:

  * `farthest_point` must SPAN a set -- it is what decides which members of a
    9,600-anchor clade step 4 ever sees.
  * `hash_sample` must be REPRESENTATIVE -- it backs a median-identity estimate,
    where a diversity-spread sample would be biased low by construction.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from protein_ncrna.diversity import (farthest_point, hash_sample, jaccard,
                                     kmers, nearest_representatives,
                                     stratified_sample)

AA = "ACDEFGHIKLMNPQRSTVWY"


def _clade(n_groups=4, per_group=25, length=200, seed=0):
    """Synthetic clade: tight subgroups around well-separated founders."""
    rnd = random.Random(seed)
    seqs = {}
    for g in range(n_groups):
        founder = "".join(rnd.choice(AA) for _ in range(length))
        for m in range(per_group):
            s = list(founder)
            for _ in range(3):                      # a few point changes only
                s[rnd.randrange(length)] = rnd.choice(AA)
            seqs[f"g{g}_m{m:03d}"] = "".join(s)
    return seqs


def test_hash_sample_is_order_independent_and_deterministic():
    ids = [f"id{i:04d}" for i in range(500)]
    a = hash_sample(ids, 50)
    b = hash_sample(list(reversed(ids)), 50)
    c = hash_sample(random.Random(7).sample(ids, len(ids)), 50)
    assert a == b == c
    assert len(a) == 50 and len(set(a)) == 50
    # A different seed gives a different sample, so the sample is not a
    # disguised alphabetical slice.
    assert hash_sample(ids, 50, seed=2) != a


def test_hash_sample_is_not_the_head_of_the_list():
    ids = [f"id{i:04d}" for i in range(500)]
    assert hash_sample(ids, 50) != sorted(ids)[:50]


def test_farthest_point_is_order_independent():
    seqs = _clade()
    ids = sorted(seqs)
    a = farthest_point(ids, seqs, 8)
    b = farthest_point(list(reversed(ids)), seqs, 8)
    c = farthest_point(random.Random(3).sample(ids, len(ids)), seqs, 8)
    assert a == b == c and len(a) == 8


def test_farthest_point_covers_every_subgroup():
    """The property the cap exists for: a slice of input order can miss whole
    subgroups, a diversity cap cannot."""
    seqs = _clade(n_groups=4, per_group=25)
    ids = sorted(seqs)
    picked = farthest_point(ids, seqs, 8)
    assert {i.split("_")[0] for i in picked} == {"g0", "g1", "g2", "g3"}
    # The old behaviour, for contrast: the first 8 in order are all one group.
    assert len({i.split("_")[0] for i in ids[:8]}) == 1


def test_farthest_point_respects_the_preference_key():
    seqs = _clade(n_groups=2, per_group=20)
    ids = sorted(seqs)
    # Mark one member of each group as complete; with everything else equal the
    # complete ones must be chosen.
    good = {ids[0], ids[25]}
    prefer = {i: (i in good,) for i in ids}
    picked = farthest_point(ids, seqs, 2, prefer=prefer)
    assert good & set(picked)


def test_cap_at_or_above_size_keeps_everything():
    seqs = _clade(n_groups=2, per_group=5)
    ids = sorted(seqs)
    assert farthest_point(ids, seqs, len(ids)) == ids
    assert farthest_point(ids, seqs, len(ids) + 10) == ids
    assert hash_sample(ids, len(ids) + 10) == ids


def test_missing_proteins_still_fill_the_cap():
    """An anchor with no protein cannot be placed in protein space, but it must
    not silently shrink the cap either."""
    seqs = _clade(n_groups=2, per_group=10)
    ids = sorted(seqs) + [f"noprot{i}" for i in range(10)]
    picked = farthest_point(ids, seqs, 12)
    assert len(picked) == 12


def test_representatives_are_not_the_outliers():
    """nearest_representatives must return typical members, not the extremes
    that seeded the groups."""
    seqs = _clade(n_groups=3, per_group=30)
    ids = sorted(seqs)
    reps = nearest_representatives(ids, seqs, 3)
    assert len(reps) == 3
    assert {i.split("_")[0] for i in reps} == {"g0", "g1", "g2"}
    for r in reps:
        mates = [i for i in ids if i.split("_")[0] == r.split("_")[0]]
        mean_to_own = sum(jaccard(kmers(seqs[r]), kmers(seqs[m]))
                          for m in mates) / len(mates)
        worst = min(mates, key=lambda m: sum(
            jaccard(kmers(seqs[m]), kmers(seqs[x])) for x in mates))
        mean_worst = sum(jaccard(kmers(seqs[worst]), kmers(seqs[m]))
                         for m in mates) / len(mates)
        assert mean_to_own >= mean_worst


def test_jaccard_and_kmers_edge_cases():
    assert jaccard(set(), {"abc"}) == 0.0
    assert jaccard({"a"}, {"a"}) == 1.0
    # A sequence shorter than k is its own token, not an empty set, so two
    # short proteins remain comparable.
    assert kmers("AC", 5) == {"AC"}
    assert len(kmers("ACDEFG", 5)) == 2


def test_stratified_sample_fills_the_cap_exactly():
    """A cap that silently returns fewer than asked is the worst of both worlds.

    The first implementation used farthest-point seeds and kept one
    representative per *non-empty* group. On near-identical proteins every
    Jaccard ties, `max()` returned the first seed every time, the groups
    collapsed, and a 648-member clade capped to 500 came back with 34 windows --
    an 93% silent data loss that looked like a working diversity cap.
    """
    seqs = _clade(n_groups=5, per_group=120)
    ids = sorted(seqs)
    for cap in (10, 50, 300, 599):
        sel = stratified_sample(ids, seqs, cap)
        assert len(sel) == cap, f"cap {cap} returned {len(sel)}"
        assert len(set(sel)) == cap


def test_stratified_sample_spreads_across_subgroups():
    seqs = _clade(n_groups=5, per_group=120)
    ids = sorted(seqs)
    sel = stratified_sample(ids, seqs, 50)
    counts = {g: sum(1 for i in sel if i.startswith(g)) for g in
              ("g0", "g1", "g2", "g3", "g4")}
    assert all(c >= 5 for c in counts.values()), counts


def test_stratified_sample_prefers_typical_members_over_outliers():
    """The distinction from farthest_point, measured rather than asserted.

    Both samplers span the set; only one of them picks *typical* members. With a
    cap well below the set size, the mean similarity of a selected member to its
    own subgroup must be higher for stratified than for farthest-point -- that
    difference is the whole reason the cap uses one and the references use the
    other.
    """
    seqs = _clade(n_groups=2, per_group=100)
    ids = sorted(seqs)

    def typicality(sel):
        tot = 0.0
        for i in sel:
            mates = [m for m in ids if m.split("_")[0] == i.split("_")[0]]
            probe = mates[::7]
            tot += sum(jaccard(kmers(seqs[i]), kmers(seqs[m]))
                       for m in probe) / len(probe)
        return tot / len(sel)

    fp, st = farthest_point(ids, seqs, 40), stratified_sample(ids, seqs, 40)
    assert typicality(st) > typicality(fp), (typicality(st), typicality(fp))
    # Both still cover both subgroups; breadth is not what differs.
    for sel in (fp, st):
        assert {i.split("_")[0] for i in sel} == {"g0", "g1"}


def test_stratified_sample_is_order_independent():
    seqs = _clade(n_groups=3, per_group=40)
    ids = sorted(seqs)
    a = stratified_sample(ids, seqs, 20)
    b = stratified_sample(list(reversed(ids)), seqs, 20)
    c = stratified_sample(random.Random(5).sample(ids, len(ids)), seqs, 20)
    assert a == b == c


def test_representatives_do_not_collapse_on_identical_sequences():
    """Exact ties must spread across seeds, not pile onto the first."""
    seqs = {f"id{i:03d}": "ACDEFGHIKLMNPQRSTVWY" * 10 for i in range(200)}
    ids = sorted(seqs)
    assert len(stratified_sample(ids, seqs, 50)) == 50
    assert len(nearest_representatives(ids, seqs, 3)) == 3


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
    print("FAILED" if fails else "all diversity tests passed")
    raise SystemExit(1 if fails else 0)
