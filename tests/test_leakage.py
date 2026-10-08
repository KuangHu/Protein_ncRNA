#!/usr/bin/env python
"""Every config must contain protein HMMs and nothing else.

The blind is the experiment. A single RNA model in an anchor or neighbourhood
config would hand the pipeline its own answer, and the ground-truth controls
would prove nothing. This checks the configs by content rather than by trusting
that whoever edited them remembered the rule.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# Directory decides the rule: configs/blind/ is loaded by blind discovery and
# must stay clean; configs/known_retron/ is the benchmark branch and may name
# msr, msd, Rfam and covariance models freely.
BLIND_DIR = ROOT / "configs/blind"
KNOWN_DIR = ROOT / "configs/known_retron"
CONFIGS = sorted(str(p.relative_to(ROOT)) for p in BLIND_DIR.glob("*.json"))

# The blind chain: every script that runs before a prediction is frozen. None of
# these may reference the benchmark branch, whatever it is named later.
BLIND_CHAIN = [
    "genomes_db_report.py", "collect_anchors.py", "refilter_census.py",
    "neighborhood_scan.py", "clade_decompose.py", "extract_windows.py",
    "permutation_null.py", "discover_blocks.py", "cap_diagnostics.py",
    "export_blocks.py", "fold_blocks.py", "locus_report.py",
    "bags_to_windows.py", "expand_by_homology.py", "bag_query_proteins.py",
    "arm2_nulls.py", "arm2_null_report.py",
]

# Boundary scripts: they compute blind columns and then *append* benchmark
# annotation to the same table. They are not in BLIND_CHAIN because they do
# legitimately read the benchmark, and they are not free of it either -- the
# invariant is that every value sourced from the benchmark lands in a column
# named `benchmark_*`, so no blind column can ever be contaminated by one and a
# reader can tell the two apart by the column name alone.
BOUNDARY_SCRIPTS = ["candidate_table.py", "triage_report.py",
                    "discovery_shortlist.py"]

# Terms that only appear when someone has started describing the RNA itself.
LEAK = re.compile(
    r"(msr|msd|omega[\s_-]?rna|tracr[\s_-]?rna|cr[\s_-]?rna|bridge[\s_-]?rna"
    r"|rfam|covariance|\bcm\b|\.cm\b)",
    re.I,
)
# Pfam accessions are PFxxxxx; Rfam is RFxxxxx and must never appear.
RFAM_ACC = re.compile(r"\bRF\d{5}\b", re.I)


def _walk(obj, path="$"):
    """Yield (json path, string) for every string in the document."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield f"{path}.{k}", str(k)
            yield from _walk(v, f"{path}.{k}")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from _walk(v, f"{path}[{i}]")
    elif isinstance(obj, str):
        yield path, obj


def test_blind_config_dir_is_clean():
    """Every file under configs/blind/ must be free of RNA vocabulary.

    This checks raw text, including prose, because a blind config has no reason
    to discuss msr/msd or Rfam at all; the benchmark branch is where that
    belongs.
    """
    offenders = []
    for p in sorted(BLIND_DIR.rglob("*")):
        if not p.is_file():
            continue
        for i, line in enumerate(p.read_text(errors="replace").splitlines(), 1):
            m = re.search(r"(\bmsr\b|\bmsd\b|retron[_\s-]?RNA|\bRfam\b|\.cm\b)",
                          line, re.I)
            if m:
                offenders.append(f"{p.relative_to(ROOT)}:{i}: {m.group(0)!r} in "
                                 f"{line.strip()[:80]}")
    assert not offenders, ("configs/blind/ must not mention msr, msd, retron_RNA, "
                           "Rfam or *.cm:\n" + "\n".join(offenders))
    print(f"ok  configs/blind/ clean ({len(list(BLIND_DIR.glob('*.json')))} files)")


def test_known_retron_branch_is_separate():
    """The benchmark branch exists and no blind code reads it."""
    assert KNOWN_DIR.is_dir(), "configs/known_retron/ must exist as the benchmark home"
    readers = []
    for py in list((ROOT / "protein_ncrna").rglob("*.py")):
        t = py.read_text()
        if "known_retron" in t or "data/known_retrons" in t:
            readers.append(str(py.relative_to(ROOT)))
    assert not readers, ("blind discovery code must not read the benchmark "
                         "branch: " + ", ".join(readers))
    # Every script in the blind chain must stay clean -- not just the four that
    # happened to be listed when this test was written. A new script that joins
    # the chain and reads the benchmark is exactly the leak this guards against,
    # and an allowlist of names would not have caught it.
    for name in BLIND_CHAIN:
        p = ROOT / "scripts" / name
        assert p.exists(), f"stale BLIND_CHAIN entry {name}"
        assert "known_retron" not in p.read_text(), \
            f"{name} is in the blind chain and reads the benchmark branch"
    # Scripts that legitimately read annotation must say so in their name or
    # docstring, so "does this file see the answers?" is answerable by reading
    # the top of it rather than by tracing imports.
    for p in sorted((ROOT / "scripts").glob("*.py")):
        t = p.read_text()
        if ("known_retron" not in t or p.name in BLIND_CHAIN
                or p.name in BOUNDARY_SCRIPTS):
            continue
        head = t[:t.find("\n\n", t.find('"""'))] if '"""' in t else ""
        assert ("UNBLIND" in head or "benchmark" in p.name
                or "known_retron" in p.name or "annotate" in p.name), (
            f"{p.name} reads the benchmark branch but is not marked UNBLIND "
            "in its docstring, named as a benchmark script, or listed in "
            "BOUNDARY_SCRIPTS")
    print(f"ok  known-retron benchmark isolated from {len(BLIND_CHAIN)} "
          "blind-chain scripts")


def test_boundary_scripts_quarantine_benchmark_columns():
    """Annotation may be appended to a blind table; it may not be mixed into it.

    `candidate_table.py` and `triage_report.py` both compute every blind column
    and then join benchmark annotation onto the same rows. That is fine, and it
    is also exactly where a leak would be invisible: one `r["score"] = ...`
    sourced from the benchmark dict and the discovery score is no longer blind.

    The invariant is structural, so it can be checked structurally: every
    assignment whose value reads the benchmark annotation must target a column
    named `benchmark_*`.
    """
    import ast as _ast
    bad = []
    for name in BOUNDARY_SCRIPTS:
        p = ROOT / "scripts" / name
        assert p.exists(), f"stale BOUNDARY_SCRIPTS entry {name}"
        tree = _ast.parse(p.read_text())
        # Names bound to the benchmark annotation, found by the file it reads.
        for node in _ast.walk(tree):
            if not isinstance(node, _ast.Assign):
                continue
            src = _ast.unparse(node.value)
            if not any(t in src for t in ("overlaps_known_retron",
                                          "members_overlapping_known",
                                          "benchmark_overlaps_known_ncrna",
                                          "benchmark_block_overlap")):
                continue
            for tgt in node.targets:
                key = (tgt.slice.value if isinstance(tgt, _ast.Subscript)
                       and isinstance(getattr(tgt, "slice", None), _ast.Constant)
                       else None)
                if isinstance(key, str) and not key.startswith("benchmark_"):
                    bad.append(f"{name}: benchmark value assigned to "
                               f"non-benchmark column {key!r}")
    assert not bad, "\n".join(bad)
    print(f"ok  {len(BOUNDARY_SCRIPTS)} boundary scripts keep benchmark values "
          "in benchmark_* columns")


def test_no_rna_terms_in_model_names_or_paths():
    """Model names and file paths must be free of RNA vocabulary.

    Prose keys (`description`, `purpose`, `_comment`) are exempt: they are where
    the rationale is written, and the census notes legitimately say things like
    'omegaRNA is encoded downstream of tnpB'. What must stay clean is anything
    the pipeline actually loads.
    """
    prose = {"description", "purpose", "_comment"}
    offenders = []
    for rel in CONFIGS:
        cfg = json.loads((ROOT / rel).read_text())
        for jpath, val in _walk(cfg):
            if any(f".{k}" in jpath for k in prose):
                continue
            m = LEAK.search(val)
            if m:
                offenders.append(f"{rel} {jpath}: {val!r} contains {m.group(0)!r}")
    assert not offenders, "RNA vocabulary in loaded config values:\n" + "\n".join(offenders)
    print(f"ok  {len(CONFIGS)} configs free of RNA terms in loaded values")


def test_every_accession_is_pfam():
    """Anchors must be protein families, never RNA families.

    Only model-declaring configs are checked. configs/blind/ also holds
    parameter configs such as step4.json, which declare no models at all; a
    blanket "must declare models" assertion would force a fake entry into them.
    """
    offenders = []
    n_model_cfgs = 0
    for rel in CONFIGS:
        cfg = json.loads((ROOT / rel).read_text())
        groups = cfg.get("families") or cfg.get("groups") or {}
        if not groups:
            continue
        n_model_cfgs += 1
        n = 0
        for fam, spec in groups.items():
            for m in spec["models"]:
                n += 1
                acc = m["acc"]
                if RFAM_ACC.match(acc) or not acc.upper().startswith("PF"):
                    offenders.append(f"{rel} {fam}/{m['name']}: accession {acc}")
        assert n, f"{rel} declares an empty model set"
    assert n_model_cfgs, "no model-declaring config found at all"
    assert not offenders, "non-Pfam accessions:\n" + "\n".join(offenders)
    print(f"ok  every accession in {n_model_cfgs} model configs is a Pfam "
          f"protein family")


def test_hmm_db_is_protein():
    n = 0
    for rel in CONFIGS:
        cfg = json.loads((ROOT / rel).read_text())
        db = cfg.get("hmm_db", "")
        if not db:
            continue  # parameter-only config, declares no database
        n += 1
        assert "Pfam" in db, f"{rel}: hmm_db {db!r} is not Pfam"
        assert not RFAM_ACC.search(db) and "Rfam" not in db, f"{rel}: hmm_db points at Rfam"
    assert n, "no config declares an HMM database"
    print(f"ok  all {n} HMM databases referenced are protein databases")


def _docstrings(tree) -> set[int]:
    """Line numbers of docstrings, which are prose and may discuss RNA tools."""
    import ast

    lines = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                             ast.AsyncFunctionDef)):
            doc = ast.get_docstring(node, clean=False)
            if doc and node.body and isinstance(node.body[0], ast.Expr):
                c = node.body[0].value
                if isinstance(c, ast.Constant) and isinstance(c.value, str):
                    lines.update(range(c.lineno, (c.end_lineno or c.lineno) + 1))
    return lines


def test_discovery_code_has_no_covariance_tooling():
    """Discovery must not invoke RNA tooling; step 6 is a separate, later stage.

    Checked against string literals in the AST rather than raw text, so that
    documenting why R-scape is deferred does not trip the test while an actual
    `which("RNAalifold")` call would.

    STEP6_SCRIPTS is an explicit allowlist, not an exemption for convenience.
    The rule being enforced is that nothing which *selects or scores* a
    candidate -- anchors, clades, windows, blocks -- may consult RNA structure
    tooling. Step 6 consumes an already-frozen candidate list and is the stage
    where folding is the whole point. Any script added here must take its input
    from a file produced upstream and must not write back into it.
    """
    import ast

    STEP6_SCRIPTS = {"fold_blocks.py"}

    rna_tools = re.compile(r"^(cmbuild|cmsearch|cmalign|cmcalibrate|R-scape|"
                           r"RNAalifold|RNAfold)$", re.I)
    offenders = []
    for py in list((ROOT / "protein_ncrna").rglob("*.py")) + \
              list((ROOT / "scripts").glob("*.py")):
        if py.name in STEP6_SCRIPTS:
            continue
        tree = ast.parse(py.read_text())
        skip = _docstrings(tree)
        for node in ast.walk(tree):
            if (isinstance(node, ast.Constant) and isinstance(node.value, str)
                    and node.lineno not in skip and rna_tools.match(node.value.strip())):
                offenders.append(f"{py.relative_to(ROOT)}:{node.lineno}: "
                                 f"{node.value!r}")
    assert not offenders, ("discovery code must not call RNA tooling:\n"
                           + "\n".join(offenders))
    # The allowlist must not silently cover a script that no longer exists.
    for name in STEP6_SCRIPTS:
        assert (ROOT / "scripts" / name).exists(), f"stale allowlist entry {name}"
    print("ok  no RNA covariance tooling invoked in discovery code")




def test_compute_scripts_require_slurm():
    """Heavy steps must refuse to run on a login node.

    Login nodes here are shared and contended; an early pilot measured ~11% CPU
    per worker, making a login-node run ~10x slower as well as antisocial.
    """
    import re as _re
    need = ["clade_decompose.py", "collect_anchors.py", "neighborhood_scan.py",
            "extract_windows.py", "refilter_census.py", "annotate_known_retrons.py",
            "genomes_db_report.py", "discover_blocks.py", "export_blocks.py",
            "locus_report.py", "triage_report.py", "cap_diagnostics.py",
            "bags_to_windows.py", "expand_by_homology.py"]
    missing = [n for n in need
               if "require_slurm" not in (ROOT / "scripts" / n).read_text()]
    assert not missing, f"compute scripts without a SLURM guard: {missing}"
    # Every such script must also have an sbatch wrapper to point people at.
    wrappers = {p.stem for p in (ROOT / "sbatch").glob("*.sbatch")}
    assert wrappers, "no sbatch wrappers found"
    print(f"ok  {len(need)} compute scripts guarded, {len(wrappers)} sbatch wrappers")


if __name__ == "__main__":
    test_blind_config_dir_is_clean()
    test_known_retron_branch_is_separate()
    test_boundary_scripts_quarantine_benchmark_columns()
    test_no_rna_terms_in_model_names_or_paths()
    test_every_accession_is_pfam()
    test_hmm_db_is_protein()
    test_discovery_code_has_no_covariance_tooling()
    test_compute_scripts_require_slurm()
    print("all leakage tests passed")
