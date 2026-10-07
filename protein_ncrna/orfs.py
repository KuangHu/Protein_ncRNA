"""Step 1a: call ORFs with prodigal and recover their genomic coordinates.

Prodigal is run in `meta` mode so that short contigs and plasmids are handled
without per-genome training. Coordinates come from the GFF (which carries the
true contig id in column 1); protein sequences come from the paired FAA. The two
files are emitted in the same order, so they are zipped rather than joined on
prodigal's internal ids, which rewrite contig names containing whitespace.
"""

from __future__ import annotations

import dataclasses
import subprocess
import tempfile
from pathlib import Path

from .seqio import read_fasta
from .tools import which


class EmptyGenome(Exception):
    """The mirrored FASTA is an empty stub rather than an assembly."""


@dataclasses.dataclass(frozen=True)
class ORF:
    contig: str
    start: int  # 1-based inclusive, forward-strand coordinates
    end: int  # 1-based inclusive
    strand: int  # +1 or -1
    index: int  # 0-based position in the prodigal output
    aa: str
    partial: str = "00"  # prodigal partial flag: '10'/'01'/'11' = runs off a contig edge

    @property
    def is_partial(self) -> bool:
        return self.partial != "00"

    @property
    def orf_id(self) -> str:
        return f"{self.contig}:{self.start}-{self.end}:{'+' if self.strand > 0 else '-'}"


def _parse_gff(path: Path) -> list[tuple[str, int, int, int, str]]:
    """Return (contig, start, end, strand, partial) per CDS.

    `partial` is prodigal's own flag for a gene running off a contig edge; such
    a gene legitimately lacks a start or stop codon, so it must be distinguished
    from a mis-oriented one rather than rediscovered downstream.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            partial = "00"
            for kv in f[8].split(";"):
                if kv.startswith("partial="):
                    partial = kv.split("=", 1)[1]
                    break
            rows.append((f[0], int(f[3]), int(f[4]), 1 if f[6] == "+" else -1, partial))
    return rows


def call_orfs(genome_fna: str | Path, workdir: str | Path | None = None) -> list[ORF]:
    """Run prodigal on one genome FASTA (may be gzipped) and return its ORFs."""
    genome_fna = Path(genome_fna)
    # genomes_db is 84% empty download stubs (20-byte gzip files). Prodigal
    # exits 18 on these, which is indistinguishable from a real parse failure,
    # so catch them here and report them as what they are.
    if genome_fna.stat().st_size < 1024:
        raise EmptyGenome(f"{genome_fna} is {genome_fna.stat().st_size} B — empty download stub")

    tmp = tempfile.TemporaryDirectory(dir=workdir) if workdir else tempfile.TemporaryDirectory()
    with tmp as td:
        td = Path(td)
        src = td / "genome.fna"
        if genome_fna.suffix == ".gz":
            with open(src, "wb") as out:
                subprocess.run(["gunzip", "-c", str(genome_fna)], stdout=out, check=True)
        else:
            src = genome_fna

        faa, gff = td / "orfs.faa", td / "orfs.gff"
        # These are single-organism assemblies, so `single` is both ~3x faster
        # than `meta` (19 s vs 53 s on a 7 Mb genome) and better suited: `meta`
        # re-scores against precomputed profiles per sequence. `single` needs
        # enough sequence to train on, so fall back to `meta` for assemblies
        # that are too short or too fragmented.
        for mode in ("single", "meta"):
            proc = subprocess.run(
                [which("prodigal"), "-i", str(src), "-a", str(faa), "-f", "gff",
                 "-o", str(gff), "-p", mode, "-q"],
                stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
            )
            if proc.returncode == 0:
                break
        else:
            raise subprocess.CalledProcessError(
                proc.returncode, proc.args,
                stderr=proc.stderr.decode(errors="replace")[-500:],
            )

        coords = _parse_gff(gff)
        prots = [s.rstrip("*") for _, s in read_fasta(str(faa))]
        if len(coords) != len(prots):
            raise RuntimeError(
                f"prodigal GFF/FAA mismatch for {genome_fna}: "
                f"{len(coords)} CDS vs {len(prots)} proteins"
            )
        # 537 assemblies in this mirror contain every contig twice (concatenated
        # downloads), so prodigal calls each gene twice and the duplicates carry
        # identical coordinates. Collapse them here: a locus present twice in a
        # file is one locus, and counting it twice inflates every downstream
        # statistic.
        orfs, seen = [], set()
        for i, ((c, s, e, st, pt), aa) in enumerate(zip(coords, prots)):
            key = (c, s, e, st)
            if key in seen:
                continue
            seen.add(key)
            orfs.append(ORF(contig=c, start=s, end=e, strand=st, index=i, aa=aa,
                            partial=pt))
        return orfs
