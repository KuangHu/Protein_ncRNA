#!/usr/bin/env python
"""Divergence and redundancy summary per anchor family.

The census question is not "how many anchors" but "how many *independent,
diverged* anchors". 13,000 IS110 hits across 13,000 near-identical K. pneumoniae
strains is one data point wearing a large hat, and step 6 would find no
covariation in it. Dereplicating at a ladder of identity thresholds separates
the two cases directly: a family whose count barely falls from 100% to 70%
identity is genuinely diverse; one that collapses is a strain artifact.

    census_report.py --census out/census --out out/census/divergence

Reads `<census>/proteins/<family>.faa` and `<census>/anchors.tsv`. Uses MMseqs2
when available; without it, writes the per-family FASTA list and leaves the
clustering columns empty rather than guessing.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_ncrna.seqio import read_fasta
from protein_ncrna.tools import which

# A family that still has many clusters at 70% identity has real phylogenetic
# spread; one that collapses to a handful is near-duplicate strains.
IDENTITY_LADDER = [1.00, 0.95, 0.90, 0.70]
# Below this many nonredundant sequences at 90% identity, a comparative RNA test
# has no power, so the family is not worth carrying into steps 3-6.
MIN_NR90_FOR_COVARIATION = 20


def mmseqs_cluster(faa: Path, ident: float, tmp: Path, mmseqs: str, threads: int) -> int | None:
    """Number of clusters at the given sequence-identity threshold."""
    pref = tmp / f"c{int(ident * 100)}"
    pref.parent.mkdir(parents=True, exist_ok=True)
    cmd = [mmseqs, "easy-cluster", str(faa), str(pref), str(tmp / f"t{int(ident*100)}"),
           "--min-seq-id", str(ident), "-c", "0.8", "--cov-mode", "0",
           "--threads", str(threads), "-v", "1"]
    proc = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        print(f"    mmseqs failed at {ident}: "
              f"{proc.stderr.decode(errors='replace')[-200:]}", file=sys.stderr)
        return None
    rep = Path(str(pref) + "_rep_seq.fasta")
    if not rep.exists():
        return None
    return sum(1 for _ in read_fasta(str(rep)))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", required=True, type=Path, help="collect_anchors --out dir")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()
    out = args.out or (args.census / "divergence")
    out.mkdir(parents=True, exist_ok=True)

    try:
        mmseqs = which("mmseqs")
    except FileNotFoundError:
        mmseqs = None
        print("mmseqs not found — clustering columns will be empty (stub).", file=sys.stderr)

    # Per-family species / genome breakdown from the anchor table.
    anchors = args.census / "anchors.tsv"
    sp_by_fam: dict[str, Counter] = defaultdict(Counter)
    genomes_by_fam: dict[str, set] = defaultdict(set)
    complete_by_fam: dict[str, list] = defaultdict(list)
    if anchors.exists():
        with open(anchors) as fh:
            cols = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                fam = f[cols["family"]]
                sp_by_fam[fam][f[cols["species"]]] += 1
                genomes_by_fam[fam].add(f[cols["genome"]])
                complete_by_fam[fam].append(f[cols["window_complete"]] == "True")

    prot_dir = args.census / "proteins"
    families = sorted(p.stem for p in prot_dir.glob("*.faa")) if prot_dir.is_dir() else []
    if not families:
        print(f"no per-family FASTA in {prot_dir}", file=sys.stderr)
        return 1

    tmp = out / "_mmseqs_tmp"
    report = []
    for fam in families:
        faa = prot_dir / f"{fam}.faa"
        n = sum(1 for _ in read_fasta(str(faa)))
        row = {
            "family": fam,
            "anchors": n,
            "genomes": len(genomes_by_fam.get(fam, ())),
            "species": len(sp_by_fam.get(fam, ())),
            "top_species": sp_by_fam[fam].most_common(3) if fam in sp_by_fam else [],
            "complete_window_frac": (
                round(sum(complete_by_fam[fam]) / len(complete_by_fam[fam]), 4)
                if complete_by_fam.get(fam) else None),
            "fasta": str(faa),
        }
        print(f"[{fam}] {n} anchors, {row['genomes']} genomes, {row['species']} species")
        for ident in IDENTITY_LADDER:
            key = f"nr{int(ident * 100)}"
            row[key] = mmseqs_cluster(faa, ident, tmp / fam, mmseqs, args.threads) \
                if mmseqs else None
            if row[key] is not None:
                print(f"    {key}: {row[key]}")
        nr90 = row.get("nr90")
        row["redundancy_90"] = round(n / nr90, 2) if nr90 else None
        row["covariation_feasible"] = (nr90 is not None and nr90 >= MIN_NR90_FOR_COVARIATION)
        report.append(row)

    if tmp.exists():
        shutil.rmtree(tmp, ignore_errors=True)

    cols = ["family", "anchors", "genomes", "species", "nr100", "nr95", "nr90", "nr70",
            "redundancy_90", "complete_window_frac", "covariation_feasible"]
    with open(out / "divergence_by_family.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(report, key=lambda r: -(r.get("nr90") or 0)):
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")
    (out / "divergence_by_family.json").write_text(json.dumps(report, indent=2))

    print(f"\n{'family':<14}{'anchors':>9}{'nr90':>8}{'redund':>8}{'species':>9}  covariation?")
    for r in sorted(report, key=lambda r: -(r.get("nr90") or 0)):
        nr90 = r.get("nr90")
        print(f"{r['family']:<14}{r['anchors']:>9}{nr90 if nr90 is not None else '-':>8}"
              f"{r.get('redundancy_90') or '-':>8}{r['species']:>9}  "
              f"{'YES' if r['covariation_feasible'] else 'no'}")
    print(f"\n(covariation? = >= {MIN_NR90_FOR_COVARIATION} nonredundant sequences at 90% id)")
    print(f"report -> {out}/divergence_by_family.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
