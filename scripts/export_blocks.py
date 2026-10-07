#!/usr/bin/env python
"""Step-5 input: core and padded FASTAs for each step-4 candidate block.

    export_blocks.py --candidates <dir>/step4_candidates.tsv \\
                     --windows <windows_dir> --out <dir>/fasta \\
                     [--soft-orf-padding-bp 75] [--min-score 0.90] [--exclude-benchmark-only]

Two files per block:

  * `<block_id>.core.fa`   -- the exact scored block, byte for byte what step 4
    called. Nothing is added, so a fold on this file is a fold on the evidence.
  * `<block_id>.padded.fa` -- the same block plus up to `--pad` bp each side.

Padding exists because a block's edges are set by the coding mask and by BLAST
HSP boundaries, neither of which is an RNA boundary. A structured RNA routinely
runs a little past where recurrence stops being detectable, and folding a
hard-trimmed block truncates the closing helix.

Padding is *soft*: it crosses predicted ORF boundaries rather than stopping at
them, and only the window edge clips it. Clipping at ORFs made the feature
vacuous -- step-4 blocks are maximal noncoding runs, so they are ORF-bounded by
construction and an ORF-clipped pad is always 0 bp. The justification for
crossing is that `prodigal -p meta` on an 11 kb fragment places a start codon
approximately, and a predicted ORF edge is not an RNA edge.

Overrun bases are lowercased and counted (`left_orf_overlap_bp`,
`right_orf_overlap_bp`, `padded_coding_frac`) so a fold that leans on them is
visible. The division of labour is strict: `padded.fa` is for exploratory
folding only, while the composite score and the permutation null are computed on
the noncoding core and are untouched by anything here.

Blind: padding uses ORF masks and window geometry only.
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_ncrna.compute import require_slurm

_d = __import__("importlib.util", fromlist=["util"])
_spec = _d.spec_from_file_location("_db", Path(__file__).resolve().parent / "discover_blocks.py")
_db = _d.module_from_spec(_spec)
_spec.loader.exec_module(_db)


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        cols = next(fh).rstrip("\n").split("\t")
        return [dict(zip(cols, line.rstrip("\n").split("\t"))) for line in fh]


def soft_pad(start: int, end: int, mask: list[tuple[int, int]], win_len: int,
             pad: int) -> dict:
    """Extend [start, end) by `pad` on each side, clipped only by the window.

    Coordinates here are *window* coordinates, not anchor-frame, because that is
    the frame the mask and the sequence live in. Predicted ORFs do not stop the
    extension; they are measured, so the caller can lowercase them.
    """
    lo = max(0, start - pad)
    hi = min(win_len, end + pad)
    l_ovl = sum(max(0, min(me, start) - max(ms, lo)) for ms, me in mask)
    r_ovl = sum(max(0, min(me, hi) - max(ms, end)) for ms, me in mask)
    return {"lo": lo, "hi": hi,
            "left_orf_overlap_bp": l_ovl, "right_orf_overlap_bp": r_ovl,
            "total_orf_overlap_bp": l_ovl + r_ovl,
            "left_window_clip_bp": pad - (start - lo),
            "right_window_clip_bp": pad - (hi - end),
            "padded_coding_frac": round((l_ovl + r_ovl) / max(1, hi - lo), 4)}


def soft_case(seq: str, lo: int, hi: int, mask: list[tuple[int, int]]) -> str:
    """Upper-case the slice, then lower-case every base inside a predicted ORF."""
    chars = list(seq[lo:hi].upper())
    for ms, me in mask:
        for i in range(max(ms, lo), min(me, hi)):
            chars[i - lo] = chars[i - lo].lower()
    return "".join(chars)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True, type=Path)
    ap.add_argument("--windows", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--soft-orf-padding-bp", "--pad", dest="pad", type=int,
                    default=75,
                    help="bp added each side of the core; crosses predicted ORFs "
                         "(overrun lowercased), clipped only at the window edge")
    ap.add_argument("--min-score", type=float, default=0.0)
    ap.add_argument("--exclude-benchmark-only", action="store_true")
    ap.add_argument("--members", type=Path,
                    help="block_members.tsv; defaults next to the candidates file")
    ap.add_argument("--tmp", type=Path, default=None)
    args = ap.parse_args()
    require_slurm("candidate block export")
    args.out.mkdir(parents=True, exist_ok=True)

    cands = [c for c in read_tsv(args.candidates)
             if float(c["score"]) >= args.min_score
             and not (args.exclude_benchmark_only and c["benchmark_only"] == "True")]
    if not cands:
        raise SystemExit("no candidates passed the filters")
    print(f"{len(cands)} candidate blocks")

    mem_path = args.members
    if mem_path is None:
        # fasta_path points into the step-4 output dir, which holds the members.
        mem_path = Path(cands[0]["fasta_path"]).parent.parent / "block_members.tsv"
    members_by_block: dict[str, list[dict]] = defaultdict(list)
    for m in read_tsv(mem_path):
        members_by_block[m["block_id"]].append(m)
    print(f"member coordinates from {mem_path}")

    by_clade: dict[str, list[dict]] = defaultdict(list)
    for c in cands:
        by_clade[c["clade_id"]].append(c)

    brows, mrows = [], []
    tmp_root = Path(tempfile.mkdtemp(dir=args.tmp))
    for clade, cs in sorted(by_clade.items()):
        recs = _db.load_clade(args.windows / f"{clade}.fasta")
        by_id = {r["anchor_id"]: r for r in recs}
        mask = _db.call_orfs_mask(recs, tmp_root / clade[:60])
        print(f"  {clade[:58]}: {len(cs)} blocks, {len(recs)} windows", flush=True)

        for c in cs:
            bid = c["block_id"]
            core_p = args.out / f"{bid}.core.fa"
            pad_p = args.out / f"{bid}.padded.fa"
            pl, pr, lo_ovl, ro_ovl, cfrac = [], [], [], [], []
            n_core = n_pad = n_clipped = 0
            with open(core_p, "w") as fc, open(pad_p, "w") as fp:
                for m in members_by_block.get(bid, []):
                    r = by_id.get(m["anchor_id"])
                    if not r:
                        continue
                    mk = mask.get(m["anchor_id"], [])
                    off = r["anchor_offset"]
                    s, e = int(m["start"]) + off, int(m["end"]) + off
                    s, e = max(0, s), min(r["win_len"], e)
                    if e - s < 10:
                        continue
                    fc.write(f">{m['anchor_id']} block={bid} "
                             f"anchor_frame={m['start']}..{m['end']}\n{r['seq'][s:e]}\n")
                    n_core += 1
                    p = soft_pad(s, e, mk, r["win_len"], args.pad)
                    lo, hi = p["lo"], p["hi"]
                    pl.append(s - lo)
                    pr.append(hi - e)
                    lo_ovl.append(p["left_orf_overlap_bp"])
                    ro_ovl.append(p["right_orf_overlap_bp"])
                    cfrac.append(p["padded_coding_frac"])
                    n_clipped += bool(p["left_window_clip_bp"]
                                      or p["right_window_clip_bp"])
                    fp.write(f">{m['anchor_id']} block={bid} "
                             f"anchor_frame={lo - off}..{hi - off} "
                             f"pad_left={s - lo} pad_right={hi - e} "
                             f"orf_overlap={p['total_orf_overlap_bp']}\n"
                             f"{soft_case(r['seq'], lo, hi, mk)}\n")
                    n_pad += 1
                    mrows.append({"block_id": bid, "anchor_id": m["anchor_id"],
                                  "core_start": m["start"], "core_end": m["end"],
                                  "pad_left_bp": s - lo, "pad_right_bp": hi - e,
                                  **{k: p[k] for k in
                                     ("left_orf_overlap_bp", "right_orf_overlap_bp",
                                      "total_orf_overlap_bp", "left_window_clip_bp",
                                      "right_window_clip_bp", "padded_coding_frac")}})
            brows.append({
                "block_id": bid, "family": c["family"], "clade_id": clade,
                "score": c["score"], "high_confidence": c["high_confidence"],
                "n_core_seqs": n_core, "n_padded_seqs": n_pad,
                "core_len": c["block_len"],
                "pad_left_bp": int(statistics.median(pl)) if pl else 0,
                "pad_right_bp": int(statistics.median(pr)) if pr else 0,
                "left_orf_overlap_bp": int(statistics.median(lo_ovl)) if lo_ovl else 0,
                "right_orf_overlap_bp": int(statistics.median(ro_ovl)) if ro_ovl else 0,
                "total_orf_overlap_bp": (int(statistics.median(lo_ovl))
                                         + int(statistics.median(ro_ovl))) if lo_ovl else 0,
                "padded_coding_frac": round(statistics.median(cfrac), 4) if cfrac else 0,
                "n_members_window_clipped": n_clipped,
                "core_fasta": str(core_p), "padded_fasta": str(pad_p),
            })

    bcols = ["block_id", "family", "clade_id", "score", "high_confidence",
             "n_core_seqs", "n_padded_seqs", "core_len", "pad_left_bp",
             "pad_right_bp", "left_orf_overlap_bp", "right_orf_overlap_bp",
             "total_orf_overlap_bp", "padded_coding_frac",
             "n_members_window_clipped", "core_fasta", "padded_fasta"]
    with open(args.out / "block_padding.tsv", "w") as fh:
        fh.write("\t".join(bcols) + "\n")
        for r in sorted(brows, key=lambda r: -float(r["score"])):
            fh.write("\t".join(str(r[c]) for c in bcols) + "\n")
    mcols = ["block_id", "anchor_id", "core_start", "core_end", "pad_left_bp",
             "pad_right_bp", "left_orf_overlap_bp", "right_orf_overlap_bp",
             "total_orf_overlap_bp", "left_window_clip_bp", "right_window_clip_bp",
             "padded_coding_frac"]
    with open(args.out / "member_padding.tsv", "w") as fh:
        fh.write("\t".join(mcols) + "\n")
        for r in mrows:
            fh.write("\t".join(str(r[c]) for c in mcols) + "\n")

    print(f"\n{len(brows)} blocks exported, {len(mrows)} member sequences")
    print(f"median pad: left {statistics.median([r['pad_left_bp'] for r in brows])}, "
          f"right {statistics.median([r['pad_right_bp'] for r in brows])} "
          f"(requested {args.pad}; only the window edge clips)")
    print(f"median ORF overrun per block: "
          f"{statistics.median([r['total_orf_overlap_bp'] for r in brows])} bp, "
          f"padded_coding_frac "
          f"{statistics.median([r['padded_coding_frac'] for r in brows]):.3f}")
    print(f"members clipped by a window edge: "
          f"{sum(r['n_members_window_clipped'] for r in brows)} / {len(mrows)}")
    print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
