"""Step 2: extract a fixed window around each anchor, oriented to anchor strand.

Orienting to the anchor strand is what makes windows from different genomes
comparable: after this step "left of the anchor" means the same thing in every
record, so a conserved block found at a consistent offset is evidence rather
than an artifact of arbitrary contig orientation. This is the same correction
the sister project had to make for insertion flanks.
"""

from __future__ import annotations

import dataclasses

from .anchors import Anchor
from .seqio import revcomp


@dataclasses.dataclass(frozen=True)
class Window:
    anchor_id: str
    genome: str
    contig: str
    contig_len: int
    family: str
    model: str
    evalue: float
    score: float
    anchor_strand: int
    # forward-strand coordinates of the extracted span, 1-based inclusive
    win_start: int
    win_end: int
    # anchor position *within the oriented window*, 0-based half-open
    anchor_offset: int
    anchor_len: int
    left_flank_bp: int
    right_flank_bp: int
    truncated_left: bool
    truncated_right: bool
    anchor_partial: bool  # anchor ORF runs off a contig edge; flanks are unreliable
    seq: str

    def to_record(self) -> dict:
        d = dataclasses.asdict(self)
        return d


def flank_availability(anchor: Anchor, contig_len: int, cfg: dict) -> dict | None:
    """Flank feasibility for one anchor, from coordinates alone.

    The census needs to know whether a usable ±flank_bp window *could* be cut,
    across tens of thousands of genomes, without paying to carve and store the
    sequence. Returns the same geometry `extract_window` would produce, oriented
    to the anchor strand, minus `seq`.
    """
    wcfg = cfg["window"]
    flank = wcfg["flank_bp"]
    if contig_len < wcfg.get("min_contig_bp", 0):
        return None

    o = anchor.orf
    win_start = max(1, o.start - flank)
    win_end = min(contig_len, o.end + flank)
    left_avail, right_avail = o.start - win_start, win_end - o.end

    if o.strand < 0 and wcfg.get("orient_to_anchor_strand", True):
        left_flank, right_flank = right_avail, left_avail
        anchor_offset = win_end - o.end
    else:
        left_flank, right_flank = left_avail, right_avail
        anchor_offset = o.start - win_start

    return {
        "win_start": win_start,
        "win_end": win_end,
        "anchor_offset": anchor_offset,
        "left_flank_bp": left_flank,
        "right_flank_bp": right_flank,
        "truncated_left": left_flank < flank,
        "truncated_right": right_flank < flank,
        "window_complete": left_flank >= flank and right_flank >= flank,
    }


def extract_window(anchor: Anchor, contig_seq: str, genome: str, cfg: dict) -> Window | None:
    """Carve ±flank_bp around one anchor and orient it to the anchor's strand.

    Returns None if the contig is shorter than the configured minimum, since a
    window that is mostly missing flank cannot support a synteny comparison.
    """
    wcfg = cfg["window"]
    flank = wcfg["flank_bp"]
    n = len(contig_seq)
    if n < wcfg.get("min_contig_bp", 0):
        return None

    o = anchor.orf
    win_start = max(1, o.start - flank)
    win_end = min(n, o.end + flank)
    seq = contig_seq[win_start - 1:win_end]

    left_avail = o.start - win_start
    right_avail = win_end - o.end

    if o.strand < 0 and wcfg.get("orient_to_anchor_strand", True):
        seq = revcomp(seq)
        anchor_offset = win_end - o.end
        left_flank, right_flank = right_avail, left_avail
    else:
        anchor_offset = o.start - win_start
        left_flank, right_flank = left_avail, right_avail

    return Window(
        anchor_id=f"{genome}|{o.orf_id}",
        genome=genome,
        contig=o.contig,
        contig_len=n,
        family=anchor.family,
        model=anchor.model,
        evalue=anchor.evalue,
        score=anchor.score,
        anchor_strand=o.strand,
        win_start=win_start,
        win_end=win_end,
        anchor_offset=anchor_offset,
        anchor_len=o.end - o.start + 1,
        left_flank_bp=left_flank,
        right_flank_bp=right_flank,
        truncated_left=left_flank < flank,
        truncated_right=right_flank < flank,
        anchor_partial=o.is_partial,
        seq=seq,
    )
