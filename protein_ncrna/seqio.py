"""FASTA reading and reverse complement. No external deps."""

from __future__ import annotations

import gzip
from typing import Iterator

_COMPLEMENT = str.maketrans("ACGTRYKMBVDHNacgtrykmbvdhn", "TGCAYRMKVBHDNtgcayrmkvbhdn")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _open(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def read_fasta(path: str) -> Iterator[tuple[str, str]]:
    """Yield (header, sequence). Header is the full line minus '>'."""
    header, chunks = None, []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            elif line:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def load_contigs(path: str) -> dict[str, str]:
    """Map first-token contig id -> sequence."""
    return {h.split()[0]: s for h, s in read_fasta(path)}
