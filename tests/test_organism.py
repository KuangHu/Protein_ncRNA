#!/usr/bin/env python
"""Organism labels must normalize to one form.

40% of real genomes in the mirror sit outside the species views and get their
label parsed from the assembly header. If that label kept its binomial casing,
'Klebsiella pneumoniae' and the view's 'klebsiella_pneumoniae' would count as
two species and inflate the diversity number the census exists to measure.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from collect_anchors import organism_from_header

CASES = [
    ("ENA|CP046278|CP046278.1 Salmonella enterica strain FDAARGOS_708 chromosome, complete genome.",
     "salmonella_enterica"),
    ("CP000026.1 Salmonella enterica subsp. enterica serovar Paratyphi A, complete genome",
     "salmonella_enterica"),
    ("NZ_CP012345.1 Klebsiella pneumoniae strain X plasmid p1", "klebsiella_pneumoniae"),
    ("JBNOAD010000088.1 Burkholderia pseudomallei isolate 42", "burkholderia_pseudomallei"),
    # Degenerate headers must not invent a species.
    ("contig_1", "unknown"),
    ("scaffold00012 length=4421 cov=12.3", "unknown"),
    (">", "unknown"),
    ("", "unknown"),
]


def test_parses_to_view_slug_form():
    for header, want in CASES:
        got = organism_from_header(header)
        assert got == want, f"{header!r} -> {got!r}, expected {want!r}"
    print(f"ok  {len(CASES)} headers parsed to view-slug form")


def test_matches_view_slugs():
    """A parsed label must be comparable to the slugs genomes_db views use."""
    view = Path("/global/scratch/users/kh36969/genomes_db/views/top10_IS110_total/by_species")
    if not view.is_dir():
        print("skip  genomes_db views not present")
        return
    slugs = {p.name for p in view.iterdir() if p.is_dir()}
    parsed = {organism_from_header(h) for h, _ in CASES}
    overlap = parsed & slugs
    assert overlap, f"no parsed label matched a view slug; parsed={parsed - {'unknown'}}"
    print(f"ok  parsed labels match real view slugs: {sorted(overlap)}")


if __name__ == "__main__":
    test_parses_to_view_slug_form()
    test_matches_view_slugs()
    print("all organism-label tests passed")
