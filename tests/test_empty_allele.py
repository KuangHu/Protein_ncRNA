"""The empty-allele rebuild turns an observation into the acceptance test.

Two things decide whether that observation means anything: what counts as an
empty allele, and where the insertion point lands after the alignment. Both are
silent when wrong -- a gapped hit still yields a sequence, and an off-by-one
insertion point still yields a 10 kb window -- so both are pinned here.
"""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def test_empty_allele_requires_one_ungapped_full_length_hsp():
    """A genome carrying the element cannot align across the junction in one
    piece. That is the whole test, so the filter must reject a gapped or
    partial HSP rather than treat it as a quality issue."""
    src = (ROOT / "scripts" / "rebuild_empty_alleles.py").read_text()
    assert "if gaps or pid < args.min_identity or qs > 3 or qe < 118:" in src
    assert "continue" in src.split(
        "if gaps or pid < args.min_identity or qs > 3 or qe < 118:")[1][:40]


def test_insertion_point_is_mapped_through_the_hsp_not_assumed():
    """insertion_point_in_flank is 60 for every record, but the HSP need not
    start at query base 1. The genomic point must be derived from the HSP
    offset, and must mirror for a minus-strand hit."""
    src = (ROOT / "scripts" / "rebuild_empty_alleles.py").read_text()
    assert "off = ip - qs" in src
    assert "gp = ss + off if plus else ss - off" in src


def test_missing_empty_allele_is_untested_not_rejected():
    """Absence of a reconstructed empty site is a fact about the genome sample,
    not about the block, so it must never become a rejection."""
    src = (ROOT / "scripts" / "empty_allele_test.py").read_text()
    assert "untested_no_empty_allele" in src
    assert "host_present_in_empty" in src
    assert "element_borne_absent_from_empty" in src


def test_both_scripts_guard_the_login_node():
    for name in ("rebuild_empty_alleles.py", "empty_allele_test.py"):
        src = (ROOT / "scripts" / name).read_text()
        assert "def require_slurm" in src, name
        assert "require_slurm()" in src.split("def main")[1], name


def test_scripts_refuse_to_run_without_slurm():
    for name in ("rebuild_empty_alleles.py", "empty_allele_test.py"):
        r = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / name),
             "--help"], capture_output=True, text=True)
        assert r.returncode == 0, f"{name} --help failed: {r.stderr}"


def test_genome_map_covers_every_species_the_targets_use():
    """A species with no genome directory yields no empty allele at all, which
    would silently shrink the tested population."""
    gm = ROOT / "data" / "genome_map.tsv"
    have = {l.split("\t")[0] for l in gm.read_text().splitlines()
            if l.strip() and not l.startswith("#")}
    need = {"ecoli", "kpneu", "efaecium", "efaecalis", "senterica",
            "paeruginosa", "sa_k24", "spneumoniae", "lmonocytogenes"}
    assert need <= have, f"missing genome sources for {need - have}"


def test_stage2_reads_side_and_strand_rather_than_assuming_them():
    """Windows are strand-normalised, so a plus hit is the tnp's direction --
    but which side the candidate sits on must be measured per window, and a
    candidate that scatters must stay visible rather than be averaged away."""
    src = (ROOT / "scripts" / "confirm_in_windows.py").read_text()
    assert 'side = "upstream" if a1 <= 0 else (' in src
    assert "frac_modal_side" in src and "frac_same_strand_as_tnp" in src


def test_stage2_null_is_same_family_other_bags():
    src = (ROOT / "scripts" / "confirm_in_windows.py").read_text()
    assert "byfam[family_of(bag)] if b != bag" in src
    assert "enrichment_over_offbag" in src


def test_stage1_emits_a_real_member_not_a_consensus():
    """Stage 2 searches real windows, so the candidate has to be a sequence
    that exists in a genome."""
    src = (ROOT / "scripts" / "nc_rna_in_insertion.py").read_text()
    assert "cand = mem[len(mem) // 2]" in src
    assert "want_members=True" in src


def test_recurrence_counts_distinct_sites_not_sequences():
    """One site contributing three near-identical regions is one observation."""
    src = (ROOT / "scripts" / "nc_rna_in_insertion.py").read_text()
    assert "k = len(set(g))" in src
    assert "rec = k / n_sites" in src
