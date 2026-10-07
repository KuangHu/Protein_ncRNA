#!/usr/bin/env python
"""Family assignment must be priority-first, not score-first.

A group II intron maturase hits GIIM and RVT_1 both. If the higher-scoring
domain won, the same biological family would scatter between GroupII_RT and RT
depending on which domain happened to score better in each sequence — and the
split would be invisible in the output.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from protein_ncrna.anchors import family_priority, load_config

CFG = load_config(ROOT / "configs/blind/anchors.json")


def _resolve(hits):
    """Mirror of find_anchors' resolution step: (family, model, evalue, score)."""
    prio = family_priority(CFG)
    fam = max(hits, key=lambda h: (prio.get(h[0], 0), h[3]))[0]
    in_fam = [h for h in hits if h[0] == fam]
    return fam, max(in_fam, key=lambda h: h[3])[1]


def test_specific_family_beats_generic_rt():
    # RVT_1 scores much higher, but GroupII_RT is the specific family.
    hits = [("RT", "RVT_1", 1e-40, 300.0), ("GroupII_RT", "GIIM", 1e-10, 60.0)]
    fam, model = _resolve(hits)
    assert fam == "GroupII_RT", f"got {fam}; generic RT won on score"
    assert model == "GIIM"
    print(f"ok  RVT_1(300) + GIIM(60) -> {fam}/{model}")


def test_score_breaks_ties_within_a_family():
    hits = [("IS110", "DEDD_Tnp_IS110", 1e-20, 50.0),
            ("IS110", "Transposase_20", 1e-30, 90.0)]
    fam, model = _resolve(hits)
    assert (fam, model) == ("IS110", "Transposase_20"), f"got {fam}/{model}"
    print(f"ok  tie within family resolved by score -> {model}")


def test_plain_rt_still_assigned_to_rt():
    hits = [("RT", "RVT_1", 1e-40, 300.0)]
    assert _resolve(hits)[0] == "RT"
    print("ok  RT-only hit stays RT")


def test_priorities_are_declared():
    prio = family_priority(CFG)
    assert prio.get("GroupII_RT", 0) > prio.get("RT", 0), \
        "GroupII_RT must outrank generic RT in configs/blind/anchors.json"
    missing = [f for f, p in prio.items() if p == 0]
    assert not missing, f"families without an explicit priority: {missing}"
    print(f"ok  priorities declared for all {len(prio)} families")




def test_models_belong_to_one_family_only():
    """model -> family is a single mapping downstream, so overlap must be refused."""
    import json, tempfile
    cfg = json.loads((ROOT / "configs/blind/anchors.json").read_text())
    fam = next(iter(cfg["families"]))
    other = [f for f in cfg["families"] if f != fam][0]
    cfg["families"][other]["models"].append(cfg["families"][fam]["models"][0])
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        json.dump(cfg, fh)
        p = fh.name
    try:
        load_config(p)
    except ValueError as e:
        assert "exactly one family" in str(e)
        print(f"ok  duplicate model rejected: {str(e)[:70]}...")
    else:
        raise AssertionError("load_config accepted a model in two families")
    finally:
        Path(p).unlink(missing_ok=True)


def test_length_floors_are_sane():
    """A family's length floor must not exceed the protein it is meant to catch."""
    f = CFG["families"]
    assert f["RT"]["min_aa_len"] == 250, "RT floor guards against the 138 aa host protein"
    assert f["TnpB_Cas12f"]["min_aa_len"] <= 400, "TnpB is ~400 aa; floor must admit it"
    assert f["Cas12"]["min_aa_len"] >= 800, "true Cas12a is ~1300 aa"
    assert "Cas12f1-like_TNB" not in [m["name"] for m in f["Cas12"]["models"]], \
        "the TNB domain is TnpB, not Cas12"
    print("ok  per-family length floors consistent with expected protein sizes")


if __name__ == "__main__":
    test_specific_family_beats_generic_rt()
    test_score_breaks_ties_within_a_family()
    test_plain_rt_still_assigned_to_rt()
    test_priorities_are_declared()
    test_models_belong_to_one_family_only()
    test_length_floors_are_sane()
    print("all family-assignment tests passed")
