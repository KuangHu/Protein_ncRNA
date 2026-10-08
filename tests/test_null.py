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
