#!/usr/bin/env python
"""Normalize the RetronDB census CSV into known_retrons.tsv.

Benchmark data only — see data/known_retrons/README.md. Nothing that produces
blind predictions may read the output.

    import_retrondb.py --csv data/known_retrons/raw/retrondb.csv \\
                       --out data/known_retrons/known_retrons.tsv

Also writes `known_retron_ncrna.fasta` next to it, holding the subset that
carries a sequence — the only subset that can actually be matched, since
RetronDB identifiers are PATRIC or protein accessions rather than assembly
accessions.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path

# The CSV ships with these spellings; keep them rather than silently guessing.
SRC = {
    "retron_name": "retron name",
    "rt_clade": "rt/clade",
    "retron_subtype": "retron (sub)",
    "msr_msd_family": "msr/msd familiy",
    "accession": "accesion",
    "species": "species/strain",
    "taxon": "taxon code",
    "ncrna": "ncRNA",
    "cluster": "cluster/domain",
}

OUT_COLS = ["retron_id", "accession", "accession_style", "retron_name", "rt_clade",
            "retron_subtype", "msr_msd_family", "cluster", "species", "taxon",
            "ncrna_len", "ncrna"]


def accession_style(acc: str) -> str:
    if not acc:
        return "none"
    if acc.startswith("fig|"):
        return "patric"
    if re.match(r"^(WP|NP|YP|AP)_\d+", acc):
        return "refseq_protein"
    if re.match(r"^GC[AF]_\d+", acc):
        return "assembly"
    if re.match(r"^[A-Z]{3}\d+\.\d+$", acc):
        return "genbank_protein"
    return "other"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    with open(args.csv, newline="", encoding="utf-8", errors="replace") as fh:
        rows = list(csv.DictReader(fh))
    missing = [v for v in SRC.values() if rows and v not in rows[0]]
    if missing:
        raise SystemExit(f"RetronDB CSV is missing expected columns: {missing}\n"
                         f"present: {list(rows[0])}")

    styles, clades, with_seq = Counter(), Counter(), 0
    out_rows = []
    for i, r in enumerate(rows):
        acc = (r.get(SRC["accession"]) or "").strip()
        nc = re.sub(r"\s+", "", (r.get(SRC["ncrna"]) or "")).upper()
        if nc and not set(nc) <= set("ACGTUN"):
            nc = ""  # not a nucleotide sequence; do not pass it to blast
        styles[accession_style(acc)] += 1
        clades[(r.get(SRC["rt_clade"]) or "").strip()] += 1
        with_seq += bool(nc)
        out_rows.append({
            "retron_id": (r.get("_id") or f"row{i}").strip(),
            "accession": acc,
            "accession_style": accession_style(acc),
            "retron_name": (r.get(SRC["retron_name"]) or "").strip(),
            "rt_clade": (r.get(SRC["rt_clade"]) or "").strip(),
            "retron_subtype": (r.get(SRC["retron_subtype"]) or "").strip(),
            "msr_msd_family": (r.get(SRC["msr_msd_family"]) or "").strip(),
            "cluster": (r.get(SRC["cluster"]) or "").strip(),
            "species": (r.get(SRC["species"]) or "").strip(),
            "taxon": (r.get(SRC["taxon"]) or "").strip(),
            "ncrna_len": len(nc),
            "ncrna": nc,
        })

    with open(args.out, "w") as fh:
        fh.write("\t".join(OUT_COLS) + "\n")
        for r in out_rows:
            fh.write("\t".join(str(r[c]) for c in OUT_COLS) + "\n")

    fasta = args.out.with_name("known_retron_ncrna.fasta")
    with open(fasta, "w") as fh:
        for r in out_rows:
            if r["ncrna"]:
                name = r["retron_name"] or r["retron_id"]
                fh.write(f">{r['retron_id']} name={name!r} clade={r['rt_clade']} "
                         f"sub={r['retron_subtype']} len={r['ncrna_len']}\n"
                         f"{r['ncrna']}\n")

    report = {
        "rows": len(out_rows),
        "with_ncrna": with_seq,
        "accession_styles": dict(styles.most_common()),
        "rt_clades": dict(clades.most_common()),
        # Stated here so no downstream recall figure can quietly use the wrong
        # denominator: only sequences can be matched, not accessions.
        "matchable_by_sequence": with_seq,
        "matchable_by_accession": 0,
        "note": ("RetronDB identifiers are PATRIC or protein accessions; this "
                 "pipeline calls ORFs de novo and genomes_db is keyed by "
                 "assembly accession, so no accession join exists. Recall must "
                 "be quoted against the %d sequence-bearing rows." % with_seq),
    }
    args.out.with_suffix(".report.json").write_text(json.dumps(report, indent=2))
    print(f"{len(out_rows)} retrons -> {args.out}")
    print(f"{with_seq} with ncRNA -> {fasta}")
    print(f"accession styles: {dict(styles.most_common())}")
    print(f"\nrecall denominator is {with_seq}, not {len(out_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
