#!/usr/bin/env python
"""The blind is the whole experiment, so it is enforced by test rather than by
discipline: no RNA model may enter an anchor config, and no discovery module may
reach into unblind/.
"""

from __future__ import annotations

import json
import re
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from protein_ncrna.anchors import load_config

# Match actual uses, not mentions: importing or opening unblind/, or hardcoding
# an Rfam accession. A word like "Rfam" inside a guard's error message is fine.
FORBIDDEN = re.compile(
    r"""(?x)
      (^|[^\w.])import\s+unblind\b
    | \bfrom\s+\.{0,2}unblind\b
    | ["'][^"']*\bunblind/
    | \bRF\d{5}\b
    """,
    re.I,
)


def test_shipped_config_is_protein_only():
    cfg = load_config(ROOT / "configs/blind/anchors.json")
    accs = [m["acc"] for spec in cfg["families"].values() for m in spec["models"]]
    assert accs, "no anchor models configured"
    assert all(a.startswith("PF") for a in accs), f"non-Pfam accession: {accs}"
    print(f"ok  {len(accs)} anchor models, all Pfam protein families")


def test_rfam_model_is_rejected():
    cfg = json.loads((ROOT / "configs/blind/anchors.json").read_text())
    cfg["families"]["retron_RNA"] = {"description": "cheating",
                                     "models": [{"name": "Retron_msr", "acc": "RF01764"}]}
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        json.dump(cfg, fh)
        p = fh.name
    try:
        load_config(p)
    except ValueError as e:
        print(f"ok  Rfam model rejected: {e}")
    else:
        raise AssertionError("load_config accepted an Rfam RNA model")
    finally:
        Path(p).unlink(missing_ok=True)


def test_discovery_code_does_not_touch_unblind():
    offenders = []
    for py in (ROOT / "protein_ncrna").rglob("*.py"):
        for i, line in enumerate(py.read_text().splitlines(), 1):
            code = line.split("#", 1)[0]
            if FORBIDDEN.search(code):
                offenders.append(f"{py.relative_to(ROOT)}:{i}: {line.strip()}")
    assert not offenders, "discovery code references blinded material:\n" + "\n".join(offenders)
    print("ok  no discovery module references unblinding material")


if __name__ == "__main__":
    test_shipped_config_is_protein_only()
    test_rfam_model_is_rejected()
    test_discovery_code_does_not_touch_unblind()
    print("all blind-guard tests passed")
