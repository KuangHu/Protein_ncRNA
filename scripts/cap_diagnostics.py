#!/usr/bin/env python
"""Per-clade diagnostics for what the window cap kept and what it cost.

    cap_diagnostics.py --windows <windows_dir> --census <corrected> \\
                       --out <windows_dir>/cap_diagnostics.tsv

A cap that spreads across diversity still removes sequences, and the quantity
that matters downstream is not how many anchors survived but how many *distinct*
proteins did. A clade that goes from 300 non-redundant proteins to 12 has been
capped into covariation-underpowered territory before step 6 ever runs, and the
only visible symptom would be R-scape reporting "underpowered" -- which looks
like biology, not sampling.

So this reports, per clade: `pre_cap_nr100` / `post_cap_nr100` and the retained
fraction, plus `post_cap_median_identity`, estimated the same way step 3 does it
-- on distinct sequences, because a median over raw anchors measures copy number
rather than divergence.

Reads `cap_manifest.tsv` written by extract_windows.py. Blind: anchors and
proteins only.
"""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_ncrna.compute import require_slurm
from protein_ncrna.diversity import hash_sample
from protein_ncrna.seqio import read_fasta
from protein_ncrna.tools import which

# Below this many distinct proteins an alignment has little chance of showing
# covariation, whatever the structure is. It is a flag, not a filter.
LOW_NR100 = 20


def median_identity(ids: list[str], seqs: dict[str, str], tmp: Path,
                    threads: int, cap: int = 300, seed: int = 1):
    """Median pairwise AA identity over distinct sequences, via mmseqs."""
    first: dict[str, str] = {}
    for m in sorted(ids):
        if m in seqs:
            first.setdefault(seqs[m], m)
    uniq = sorted(first.values())
    if len(uniq) < 2:
        return 1.0 if uniq else None
    sample = hash_sample(uniq, cap, seed)
    tmp.mkdir(parents=True, exist_ok=True)
    q = tmp / "c.faa"
    with open(q, "w") as fh:
        for m in sample:
            fh.write(f">{m}\n{seqs[m]}\n")
    res = tmp / "h.m8"
    proc = subprocess.run(
        [which("mmseqs"), "easy-search", str(q), str(q), str(res), str(tmp / "t"),
         "--threads", str(threads), "-v", "1", "--max-seqs", "400",
         "--format-output", "query,target,fident"],
        stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    if proc.returncode != 0 or not res.exists():
        return None
    vals = [float(f[2]) for f in
            (l.rstrip("\n").split("\t") for l in open(res))
            if len(f) >= 3 and f[0] != f[1]]
    return round(statistics.median(vals), 4) if vals else None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", required=True, type=Path)
    ap.add_argument("--census", required=True, type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--identity", action="store_true", default=True)
    ap.add_argument("--no-identity", dest="identity", action="store_false")
    args = ap.parse_args()
    require_slurm("cap diagnostics")
    out = args.out or args.windows / "cap_diagnostics.tsv"

    man = args.windows / "cap_manifest.tsv"
    if not man.exists():
        raise SystemExit(f"no {man}; run extract_windows.py with --max-per-clade")
    pre: dict[str, list[str]] = defaultdict(list)
    post: dict[str, list[str]] = defaultdict(list)
    reason: dict[str, str] = {}
    with open(man) as fh:
        next(fh)
        for line in fh:
            clade, aid, sel, why = line.rstrip("\n").split("\t")
            pre[clade].append(aid)
            reason[clade] = why
            if sel == "True":
                post[clade].append(aid)

    # The two window producers lay the anchor proteins out differently:
    # extract_windows.py writes one faa per *family*, bags_to_windows.py one per
    # *clade*. Read both, keyed on whichever file exists, so the distinct-protein
    # count is available on either corpus.
    seqs: dict[str, str] = {}
    names = {c.split("__")[0] for c in pre} | set(pre)
    for name in sorted(names):
        faa = args.census / "proteins" / f"{name}.faa"
        if faa.exists():
            for h, s in read_fasta(str(faa)):
                seqs.setdefault(h.split()[0], s)
    print(f"{len(pre)} clades, {len(seqs)} anchor proteins loaded")
    # Every nr100 column is derived from `seqs`. With none loaded the tool still
    # writes a full table, of zeros, and a cap that discarded 95% of a clade's
    # distinct proteins is indistinguishable from one that discarded none. That
    # is a silent wrong answer to the only question this script exists to ask.
    if not seqs:
        raise SystemExit(
            f"no anchor proteins loaded from {args.census / 'proteins'}; "
            f"expected <family>.faa or <clade>.faa for families "
            f"{sorted({c.split('__')[0] for c in pre})[:5]}. Refusing to report "
            f"nr100 counts of 0 as if the cap cost nothing.")

    nr = lambda ids: len({seqs[i] for i in ids if i in seqs})
    tmp_root = Path(tempfile.mkdtemp())
    rows = []
    for clade in sorted(pre):
        a, b = nr(pre[clade]), nr(post[clade])
        mi = (median_identity(post[clade], seqs, tmp_root / clade[:60],
                              args.threads) if args.identity else "")
        rows.append({
            "clade_id": clade, "cap_reason": reason[clade],
            "pre_cap_anchors": len(pre[clade]), "post_cap_anchors": len(post[clade]),
            "pre_cap_nr100": a, "post_cap_nr100": b,
            "nr100_retained_frac": round(b / a, 4) if a else "",
            "post_cap_median_identity": "" if mi is None else mi,
            # The two ways a clade arrives at step 6 with no covariation power:
            # too few distinct sequences, or distinct sequences that are nearly
            # the same. Both are sampling-visible before any folding happens.
            "low_power_flag": ("few_distinct" if b < LOW_NR100 else
                               "near_identical" if (mi or 0) >= 0.95 else ""),
        })
        print(f"  {clade[:54]:<54} nr100 {a} -> {b}"
              + (f", id {mi}" if mi is not None else ""), flush=True)

    cols = ["clade_id", "cap_reason", "pre_cap_anchors", "post_cap_anchors",
            "pre_cap_nr100", "post_cap_nr100", "nr100_retained_frac",
            "post_cap_median_identity", "low_power_flag"]
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    capped = [r for r in rows if r["post_cap_anchors"] < r["pre_cap_anchors"]]
    summary = {
        "clades": len(rows), "clades_capped": len(capped),
        "median_nr100_retained_frac_capped": (statistics.median(
            [r["nr100_retained_frac"] for r in capped]) if capped else None),
        "worst_nr100_loss": sorted(
            ({"clade_id": r["clade_id"], "pre": r["pre_cap_nr100"],
              "post": r["post_cap_nr100"]} for r in capped),
            key=lambda r: r["post"] / max(1, r["pre"]))[:5],
        "clades_few_distinct": sum(1 for r in rows
                                   if r["low_power_flag"] == "few_distinct"),
        "clades_near_identical": sum(1 for r in rows
                                     if r["low_power_flag"] == "near_identical"),
        "low_nr100_threshold": LOW_NR100,
    }
    Path(str(out).replace(".tsv", ".json")).write_text(json.dumps(summary, indent=2))
    print("\n" + json.dumps(summary, indent=2))
    print(f"-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
