#!/usr/bin/env python
"""Regression: clade identity must measure divergence, not copy number.

This is the bug that invalidated an entire published result set, so it gets a
behavioural test rather than a comment. Step 3 estimated each clade's median
pairwise AA identity by sampling its *raw members*. In these families a clade is
mostly exact duplicates -- one real GroupII_RT clade has 799 members and 29
distinct proteins, a ratio of 0.036 -- so the sample was overwhelmingly
duplicate-against-duplicate and the median converged on 1.0 no matter how
diverged the clade actually was.

That number is not cosmetic: it gates clade selection. The same clade measured
0.451 from an ordered sample and 1.000 from a random one, and the 1.000 crossed
the `identity > 0.9` rejection gate, which dropped GroupII_RT from the 0.35 rung
to the 0.70 rung and 3,415 alignable windows to 217. Neither number was a
measurement of anything.

The fixture reproduces that exact shape (799 raw / 29 distinct) and asserts on
what is actually handed to the aligner, because that is the only place the bug
is visible: a refactor that reintroduces raw-member sampling would still return
a plausible-looking float.
"""

from __future__ import annotations

import importlib.util
import random
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

_spec = importlib.util.spec_from_file_location(
    "clade_decompose", ROOT / "scripts" / "clade_decompose.py")
cd = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cd)

AA = "ACDEFGHIKLMNPQRSTVWY"
N_DISTINCT = 29
N_MEMBERS = 799


def _duplicate_heavy_clade(seed: int = 0):
    """799 members over 29 distinct sequences, the real clade's ratio."""
    rnd = random.Random(seed)
    distinct = ["".join(rnd.choice(AA) for _ in range(300))
                for _ in range(N_DISTINCT)]
    seqs, members = {}, []
    for i in range(N_MEMBERS):
        mid = f"anchor{i:04d}"
        members.append(mid)
        seqs[mid] = distinct[i % N_DISTINCT]      # every distinct seq recurs
    return members, seqs, distinct


class _Capture:
    """Stands in for mmseqs: records the query FASTA, returns uniform hits."""

    def __init__(self):
        self.n_query_records = None
        self.query_seqs = None

    def run(self, cmd, **kw):
        q = Path(cmd[2])
        res = Path(cmd[4])
        recs = [l for l in open(q) if l.startswith(">")]
        self.n_query_records = len(recs)
        self.query_seqs = {l.strip() for l in open(q) if not l.startswith(">")}
        ids = [l[1:].strip() for l in recs]
        with open(res, "w") as fh:
            for a in ids:
                for b in ids:
                    fh.write(f"{a}\t{b}\t0.5000\n")

        class R:
            returncode = 0
            stderr = b""
        return R()


def _call(members, seqs, cap=300):
    cap_obj = _Capture()
    real_run, real_which = cd.subprocess.run, cd.which
    cd.subprocess.run = cap_obj.run
    cd.which = lambda _n: "mmseqs"
    try:
        with tempfile.TemporaryDirectory() as td:
            val = cd.approx_median_identity(
                Path(td) / "x.faa", members, seqs, Path(td) / "t", 1, cap=cap)
    finally:
        cd.subprocess.run, cd.which = real_run, real_which
    return val, cap_obj


def test_identity_samples_distinct_sequences_not_raw_members():
    """The whole bug in one assertion: 29 go to the aligner, not 799."""
    members, seqs, distinct = _duplicate_heavy_clade()
    val, cap = _call(members, seqs)
    assert cap.n_query_records == N_DISTINCT, (
        f"identity estimate sent {cap.n_query_records} sequences to the aligner; "
        f"the clade has {N_MEMBERS} members but only {N_DISTINCT} distinct "
        "proteins, so anything above 29 is sampling copy number")
    assert len(cap.query_seqs) == N_DISTINCT
    assert val == 0.5


def test_identity_does_not_depend_on_copy_number():
    """Same 29 proteins, wildly different copy numbers -> same estimate.

    If duplicates leak back in, the two clades diverge: the one with 40x copies
    is almost all self-pairs and drifts toward 1.0.
    """
    rnd = random.Random(1)
    distinct = ["".join(rnd.choice(AA) for _ in range(300))
                for _ in range(N_DISTINCT)]
    out = []
    for copies in (1, 40):
        seqs, members = {}, []
        for i in range(N_DISTINCT * copies):
            mid = f"a{i:05d}"
            members.append(mid)
            seqs[mid] = distinct[i % N_DISTINCT]
        _, cap = _call(members, seqs)
        out.append(cap.n_query_records)
    assert out[0] == out[1] == N_DISTINCT, out


def test_single_distinct_sequence_is_perfect_identity_without_aligning():
    """A clade of one distinct protein is 1.0 by definition, not by measurement.

    It must not reach the aligner at all -- an all-vs-all of identical sequences
    is both wasteful and, with no off-diagonal pairs, undefined.
    """
    seq = "".join(random.Random(2).choice(AA) for _ in range(300))
    seqs = {f"m{i}": seq for i in range(500)}
    val, cap = _call(sorted(seqs), seqs)
    assert val == 1.0
    assert cap.n_query_records is None, "should not have invoked the aligner"


def test_identity_cap_applies_after_deduplication():
    """The cap bounds distinct sequences, so it cannot be spent on duplicates."""
    rnd = random.Random(3)
    distinct = ["".join(rnd.choice(AA) for _ in range(200)) for _ in range(500)]
    seqs, members = {}, []
    for i in range(5000):
        mid = f"m{i:05d}"
        members.append(mid)
        seqs[mid] = distinct[i % 500]
    _, cap = _call(members, seqs, cap=300)
    assert cap.n_query_records == 300
    assert len(cap.query_seqs) == 300, "the capped sample must still be distinct"


def test_fewer_than_two_members_returns_none():
    assert cd.approx_median_identity(Path("x"), ["only"], {"only": "ACDE"},
                                     Path("t"), 1) is None


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
    print("FAILED" if fails else "all identity tests passed")
    raise SystemExit(1 if fails else 0)
