#!/usr/bin/env python
"""Step 4 — recurrent noncoding blocks inside a clade, in anchor coordinates.

Blind. Nothing here reads RetronDB, Rfam, a covariance model, or any known
ncRNA family: a block is called from recurrence and geometry alone. The
benchmark comparison happens afterwards, in a different script, against output
this one produced without knowing the answer.

    discover_blocks.py --windows <dir> --out <dir> [--decoys <dir>]

Coordinate frame
----------------
`extract_windows.py` already reverse-complemented minus-strand loci, so every
window reads in the anchor's direction. This script shifts to the **anchor
frame**: x = position_in_window - anchor_offset, so x=0 is the anchor's first
base and the anchor occupies [0, anchor_len). Upstream is negative, downstream
is >= anchor_len, and the frame is comparable across members whose anchors
differ in length. Every coordinate emitted is in this frame.

Method
------
1. Prodigal (-p meta, windows are fragments) on every window; all ORFs plus the
   anchor itself become a coding mask. Blocks are called only outside it.
2. Up to K reference windows per clade. Members are blastn'd against them with
   `-task blastn -word_size 11` and `-dust no` -- megablast misses diverged
   short structured RNAs, and dust would mask exactly the AT-rich low-complexity
   stretches a real ncRNA may contain.
3. Each member keeps its best-scoring reference, so a member contributes once.
   HSPs more than half-covered by coding sequence on either side are dropped.
4. Per-base recurrence is accumulated in the anchor frame; maximal runs at or
   above `--min-recurrence` and at least `--min-len` long become blocks.
5. Each block is scored on recurrence, coordinate tightness, strand consistency
   and length consistency, and the per-member intervals are written out so a
   block can be traced back to a locus.

A high score is a candidate, not a finding. Step 6 (covariation) and step 7
(coding/repeat/UTR exclusion) are what turn one into evidence.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import statistics
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_ncrna.compute import require_slurm
from protein_ncrna.diversity import jaccard as _jaccard, kmers as _kmers, \
    nearest_representatives
from protein_ncrna.seqio import read_fasta
from protein_ncrna.tools import which

BLAST_FMT = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"


# --------------------------------------------------------------------------- #
# windows

def parse_header(h: str) -> dict:
    f = h.split()
    d = {"anchor_id": f[0]}
    for kv in f[1:]:
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    d["anchor_offset"] = int(d.get("anchor_offset", 0))
    d["anchor_len"] = int(d.get("anchor_len", 0))
    return d


def load_clade(path: Path) -> list[dict]:
    out = []
    for h, s in read_fasta(str(path)):
        d = parse_header(h)
        d["seq"] = s
        d["win_len"] = len(s)
        out.append(d)
    return out


# --------------------------------------------------------------------------- #
# coding mask

def call_orfs_mask(recs: list[dict], tmp: Path) -> dict[str, list[tuple[int, int]]]:
    """Window-coordinate coding intervals per anchor, 0-based half-open.

    Prodigal is run once over the whole clade as a multi-FASTA in meta mode;
    windows are 11 kb fragments, so single-genome training is neither available
    nor appropriate. The anchor's own interval is added unconditionally -- it is
    coding by construction and must never be called as a block.
    """
    tmp.mkdir(parents=True, exist_ok=True)
    fna = tmp / "win.fna"
    gff = tmp / "win.gff"
    # Prodigal truncates ids at whitespace; anchor ids contain '|' and ':' but
    # no spaces, so the first token round-trips.
    with open(fna, "w") as fh:
        for r in recs:
            fh.write(f">{r['anchor_id']}\n{r['seq']}\n")
    proc = subprocess.run(
        [which("prodigal"), "-i", str(fna), "-f", "gff", "-o", str(gff),
         "-p", "meta", "-q"],
        stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    mask: dict[str, list[tuple[int, int]]] = defaultdict(list)
    if proc.returncode == 0 and gff.exists():
        with open(gff) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 8 or f[2] != "CDS":
                    continue
                mask[f[0]].append((int(f[3]) - 1, int(f[4])))
    else:
        print(f"  prodigal failed: {proc.stderr.decode(errors='replace')[-160:]}",
              file=sys.stderr)
    for r in recs:
        a = r["anchor_offset"]
        mask[r["anchor_id"]].append((a, a + r["anchor_len"]))
    return {k: _merge(v) for k, v in mask.items()}


def _merge(iv: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[list[int]] = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def trim_member(a: int, b: int, ba: int, bb: int, ma: int, mb: int,
                strand: str) -> tuple[int, int] | None:
    """Cut a member interval down to the part answering to a block.

    `[a, b)` is a segment on the reference axis and `[ma, mb)` is the same
    piece of sequence on the member's axis; `[ba, bb)` is the called block.
    Both intervals are stored ascending, but a minus-strand HSP maps them in
    opposite directions, so trimming the reference's *left* end must shorten
    the member's *right* end.

    Applying the plus-strand formula to a minus HSP returns a mirrored slice --
    the right length, at the wrong place inside the HSP -- and does so
    silently, propagating into block_members.tsv, the block FASTAs, the padded
    exports and every triage coordinate. `tests/test_blocks.py` pins the
    mirror property rather than leaving it to inspection of BLAST output.

    THE TRIM IS NOT 1:1. A BLAST HSP is gapped, so the reference span and the
    member span are different lengths, and a trim measured on the reference
    axis over-cuts the member axis in proportion to the indels between them.
    Where the overlap is small and the member carries a deletion, the two
    trims together exceed the member's own span and the coordinates *cross*:
    `end < start`. That is not a rounding error, it is a member interval that
    no longer denotes any sequence, and it propagated into `median_len` (one
    block came out at -67 bp), through `length_consistency = 1 - mad/med_len`,
    which is only bounded in [0, 1] while `med_len > 0`, and into a composite
    score of 2.15 on a scale that ends at 1.

    So the trims are scaled by the HSP's own length ratio, clamped inside the
    member's HSP, and a member trimmed to nothing returns None -- it does not
    overlap the block in any meaningful sense and must be dropped rather than
    carried with crossed coordinates.
    """
    ref_span, mem_span = b - a, mb - ma
    if ref_span <= 0 or mem_span <= 0:
        return None
    left_trim = max(a, ba) - a
    right_trim = b - min(b, bb)
    scale = mem_span / ref_span
    lt, rt = int(round(left_trim * scale)), int(round(right_trim * scale))
    if strand == "+":
        fa, fb = ma + lt, mb - rt
    else:
        fa, fb = ma + rt, mb - lt
    # The trimmed interval can only ever be a sub-interval of the HSP it came
    # from; anything outside that is an artifact of the scaling.
    fa = max(ma, min(fa, mb))
    fb = max(ma, min(fb, mb))
    return (fa, fb) if fb > fa else None


def subtract(s: int, e: int, mask: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """The parts of [s, e) left after removing coding sequence.

    Trimming rather than discarding is essential: these windows are ~86% coding,
    so almost every BLAST HSP spans a gene *and* the intergenic stretch beside
    it. Dropping such an HSP whole discards the only part that could be an ncRNA
    -- in testing that is what reduced a 493-member clade to 154 usable members
    and made the recurrence threshold unreachable.
    """
    out, cur = [], s
    for ms, me in mask:
        if me <= cur:
            continue
        if ms >= e:
            break
        if ms > cur:
            out.append((cur, min(ms, e)))
        cur = max(cur, me)
        if cur >= e:
            break
    if cur < e:
        out.append((cur, e))
    return [(a, b) for a, b in out if b > a]


def coding_frac(mask: list[tuple[int, int]], s: int, e: int) -> float:
    if e <= s:
        return 1.0
    cov = sum(max(0, min(e, me) - max(s, ms)) for ms, me in mask)
    return cov / (e - s)


# --------------------------------------------------------------------------- #
# blast

_PROT: dict[str, str] = {}


def pick_references(recs: list[dict], k: int, seed: int = 1) -> list[dict]:
    """Up to k references spread across the clade's own protein diversity.

    A single reference, or three near-identical ones, leaves the diverged half
    of a clade with no BLAST hit at all -- which is silent, because those members
    simply never contribute to recurrence and the clade looks empty rather than
    under-sampled. Spreading references across protein space gives every
    subgroup something close enough to align to.

    Diversity is measured on the anchor proteins (5-mer Jaccard), never on RNA
    and never on known-RNA coordinates. Seeds are chosen farthest-point, then
    each member is assigned to its nearest seed and the group's *representative*
    is used as the reference -- the longest window, breaking ties by mean
    similarity to its own group. Taking the farthest-point seeds themselves
    would make references out of the clade's oddest members, which is the
    opposite of what a reference should be.

    Every tie is broken by a seeded hash of the anchor id rather than by list
    position, so the reference set is a property of the clade and not of the
    order its windows happened to be read in.
    """
    full = max((r["win_len"] for r in recs), default=0)
    by_id = {r["anchor_id"]: r for r in recs}
    ids = sorted(by_id)
    if len(recs) <= k:
        return sorted(recs, key=lambda r: (-r["win_len"], r["anchor_id"]))

    # Fall back to the window's own 5' sequence only when the anchor protein is
    # missing; that is a degraded comparison, not an equivalent one.
    seqs = {i: _PROT.get(i, by_id[i]["seq"][:3000]) for i in ids}
    prefer = {i: (by_id[i]["win_len"] >= full, by_id[i]["win_len"]) for i in ids}
    chosen = nearest_representatives(ids, seqs, k, prefer=prefer, seed=seed)
    return [by_id[i] for i in chosen]


def blast_clade(recs: list[dict], refs: list[dict], tmp: Path,
                threads: int, evalue: float) -> list[list[str]]:
    tmp.mkdir(parents=True, exist_ok=True)
    ref_ids = {r["anchor_id"] for r in refs}
    sub = tmp / "ref.fna"
    qry = tmp / "qry.fna"
    with open(sub, "w") as fh:
        for r in refs:
            fh.write(f">{r['anchor_id']}\n{r['seq']}\n")
    with open(qry, "w") as fh:
        for r in recs:
            if r["anchor_id"] not in ref_ids:
                fh.write(f">{r['anchor_id']}\n{r['seq']}\n")
    if qry.stat().st_size == 0:
        return []
    db = tmp / "refdb"
    subprocess.run([which("makeblastdb"), "-in", str(sub), "-dbtype", "nucl",
                    "-out", str(db)],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    res = tmp / "hits.m8"
    proc = subprocess.run(
        [which("blastn"), "-query", str(qry), "-db", str(db), "-outfmt", BLAST_FMT,
         "-evalue", str(evalue), "-num_threads", str(threads),
         "-task", "blastn", "-word_size", "11", "-dust", "no",
         "-max_target_seqs", str(max(10, len(refs) * 50))],
        stdout=open(res, "w"), stderr=subprocess.PIPE)
    if proc.returncode != 0:
        print(f"  blastn failed: {proc.stderr.decode(errors='replace')[-160:]}",
              file=sys.stderr)
        return []
    with open(res) as fh:
        return [l.rstrip("\n").split("\t") for l in fh if l.strip()]


# --------------------------------------------------------------------------- #
# blocks

def call_blocks(clade: str, recs: list[dict], args) -> tuple[list[dict], dict]:
    by_id = {r["anchor_id"]: r for r in recs}
    n = len(recs)
    stats = {"clade": clade, "n_windows": n, "n_blocks": 0}
    if n < args.min_members:
        stats["skipped"] = f"only {n} windows"
        return [], stats

    with tempfile.TemporaryDirectory(dir=args.tmp) as td:
        tmp = Path(td)
        mask = call_orfs_mask(recs, tmp / "orf")
        refs = pick_references(recs, args.references)
        rows = blast_clade(recs, refs, tmp / "blast", args.threads, args.evalue)

    ref_by_id = {r["anchor_id"]: r for r in refs}
    # Keep each member's single best reference so one member counts once.
    best_ref: dict[str, tuple[float, str]] = {}
    for f in rows:
        q, s, bits = f[0], f[1], float(f[9])
        if q not in by_id or s not in ref_by_id:
            continue
        if q not in best_ref or bits > best_ref[q][0]:
            best_ref[q] = (bits, s)

    # member -> list of (ref_frame_start, ref_frame_end, mem_frame_start,
    #                    mem_frame_end, strand, pident)
    segs: dict[str, list[tuple]] = defaultdict(list)
    for f in rows:
        q, s = f[0], f[1]
        if q not in best_ref or best_ref[q][1] != s:
            continue
        qs, qe = int(f[4]) - 1, int(f[5])
        s1, s2 = int(f[6]), int(f[7])
        # blastn reports subject coords descending on the minus strand; the
        # query side is always ascending.
        strand = "+" if s2 > s1 else "-"
        rs, re = min(s1, s2) - 1, max(s1, s2)
        qa, ra = by_id[q], ref_by_id[s]
        span_r, span_q = re - rs, qe - qs
        if span_r <= 0 or span_q <= 0:
            continue

        def to_query(x: int) -> int:
            """Position in the HSP mapped from reference to member coordinates.

            Linear interpolation across the HSP. BLAST HSPs contain gaps, so this
            is approximate; at the tens-of-bp resolution a block is called at,
            the error is small, and every coordinate is re-derived from the
            member's own sequence when the block FASTA is written.
            """
            t = (x - rs) / span_r
            return int(round(qs + t * span_q if strand == "+"
                             else qe - t * span_q))

        # A block must be noncoding in the reference *and* in the member it is
        # counted from, so trim against both masks and keep the intersection.
        for a, b in subtract(rs, re, mask.get(s, [])):
            qa0, qa1 = sorted((to_query(a), to_query(b)))
            for c, d in subtract(qa0, qa1, mask.get(q, [])):
                if d - c < args.min_len:
                    continue
                # Carry the member-side trim back onto the reference axis so the
                # two sides describe the same piece of sequence.
                ra0 = a + int((c - qa0) * span_r / max(1, qa1 - qa0))
                ra1 = a + int((d - qa0) * span_r / max(1, qa1 - qa0))
                if strand == "-":
                    ra0, ra1 = a + (b - a) - (ra1 - a), a + (b - a) - (ra0 - a)
                ra0, ra1 = max(a, min(ra0, b)), max(a, min(ra1, b))
                if ra1 - ra0 < args.min_len:
                    continue
                segs[q].append((ra0 - ra["anchor_offset"], ra1 - ra["anchor_offset"],
                                c - qa["anchor_offset"], d - qa["anchor_offset"],
                                strand, float(f[2])))

    # Reference coverage: how much of the clade the reference set actually
    # reaches. A clade with no blocks and a low recovery fraction is
    # under-referenced; one with high recovery and no blocks is a real negative.
    stats["n_reference_windows"] = len(refs)
    stats["n_members_with_any_hsp"] = len(best_ref)
    stats["frac_members_with_any_hsp"] = round(len(best_ref) / n, 4)
    stats["n_members_with_noncoding_hsp"] = len(segs)
    stats["frac_members_with_noncoding_hsp"] = round(len(segs) / n, 4)
    per_ref = Counter(best_ref[q][1] for q in segs if q in best_ref)
    stats["per_reference_recovered_members"] = {
        r["anchor_id"]: per_ref.get(r["anchor_id"], 0) for r in refs}
    stats["references_used"] = sum(1 for v in per_ref.values() if v)

    if not segs:
        stats["note"] = "no noncoding segments survived"
        return [], stats

    # Per-base recurrence in the anchor frame, counting each member once.
    cov: dict[int, set] = defaultdict(set)
    for q, ss in segs.items():
        for a, b, *_ in ss:
            for x in range(a, b):
                cov[x].add(q)
    if not cov:
        return [], stats

    need = max(args.min_members, int(round(args.min_recurrence * n)))
    lo, hi = min(cov), max(cov)
    blocks_iv, run = [], None
    for x in range(lo, hi + 1):
        if len(cov.get(x, ())) >= need:
            run = (x, x + 1) if run is None else (run[0], x + 1)
        elif run is not None:
            blocks_iv.append(run)
            run = None
    if run is not None:
        blocks_iv.append(run)
    blocks_iv = [(a, b) for a, b in blocks_iv if b - a >= args.min_len]

    out = []
    for i, (ba, bb) in enumerate(sorted(blocks_iv, key=lambda t: t[0])):
        members = []
        for q, ss in segs.items():
            best = None
            for a, b, ma, mb, strand, pid in ss:
                ov = min(b, bb) - max(a, ba)
                if ov <= 0:
                    continue
                t = trim_member(a, b, ba, bb, ma, mb, strand)
                if t is None:          # trimmed to nothing: not a member
                    continue
                fa, fb = t
                if best is None or ov > best[0]:
                    best = (ov, fa, fb, strand, pid)
            if best:
                members.append({"anchor_id": q, "start": best[1], "end": best[2],
                                "strand": best[3], "pident": round(best[4], 2),
                                "len": best[2] - best[1]})
        if len(members) < need:
            continue
        out.append(score_block(clade, i, ba, bb, members, n, by_id, best_ref))
    stats["n_blocks"] = len(out)
    if out:
        top = max(out, key=lambda b: b["composite_score"])
        stats["top_block_id"] = top["block_id"]
        stats["top_block_member_fraction"] = top["recurrence"]
    return out, stats


def _mad(xs: list[float]) -> float:
    if len(xs) < 2:
        return 0.0
    m = statistics.median(xs)
    return statistics.median([abs(x - m) for x in xs])


def score_block(clade: str, idx: int, ba: int, bb: int, members: list[dict],
                n_total: int, by_id: dict, best_ref: dict) -> dict:
    starts = [m["start"] for m in members]
    lens = [m["len"] for m in members]
    strands = [m["strand"] for m in members]

    recurrence = len(members) / n_total
    # Tightness: how reproducibly the block sits at the same distance from the
    # anchor. A real cis element holds position; a dispersed repeat does not.
    mad_start = _mad(starts)
    coord_tightness = 1.0 / (1.0 + mad_start / 50.0)
    # Windows are already anchor-oriented, so a cis-acting RNA should land on
    # one strand across the clade. A mixed-strand block is a repeat signature.
    strand_consistency = max(strands.count("+"), strands.count("-")) / len(strands)
    med_len = statistics.median(lens)
    length_consistency = max(0.0, 1.0 - (_mad(lens) / med_len)) if med_len else 0.0

    composite = (recurrence * coord_tightness * strand_consistency
                 * length_consistency)
    med_start = statistics.median(starts)
    anchor_len = statistics.median([by_id[m["anchor_id"]]["anchor_len"]
                                    for m in members])
    # A block is assembled from members that each aligned to their own best
    # reference, so name the reference that carried most of them.
    refs_used = Counter(best_ref[m["anchor_id"]][1] for m in members
                        if m["anchor_id"] in best_ref)
    return {
        "clade": clade,
        "block_id": f"{clade}__B{idx:03d}",
        "reference_id": refs_used.most_common(1)[0][0] if refs_used else "",
        "n_references_contributing": len(refs_used),
        "member_fraction": round(len(members) / n_total, 4),
        "ref_start": ba, "ref_end": bb, "ref_len": bb - ba,
        "n_members": len(members), "n_clade_windows": n_total,
        "recurrence": round(recurrence, 4),
        "coord_tightness": round(coord_tightness, 4),
        "mad_start_bp": round(mad_start, 1),
        "strand_consistency": round(strand_consistency, 4),
        "length_consistency": round(length_consistency, 4),
        "composite_score": round(composite, 4),
        "median_start": int(med_start),
        "median_len": int(med_len),
        "median_pident": round(statistics.median([m["pident"] for m in members]), 2),
        # Signed distance from the anchor: negative is upstream of the start
        # codon, positive past the stop. Reported, never used to filter.
        "dist_to_anchor": (int(med_start) if med_start < 0
                           else max(0, int(med_start - anchor_len))),
        "side": "upstream" if med_start < 0 else "downstream",
        "members": members,
    }


# --------------------------------------------------------------------------- #

def _run(task):
    path, is_decoy, args = task
    clade = path.stem
    try:
        recs = load_clade(path)
        blocks, stats = call_blocks(clade, recs, args)
        stats["decoy"] = is_decoy
        for b in blocks:
            b["decoy"] = is_decoy
            b["family"] = clade.replace("DECOY_", "").split("__")[0]
        return blocks, stats, None
    except Exception as e:
        return [], {"clade": clade, "decoy": is_decoy}, f"{type(e).__name__}: {e}"


BLOCK_COLS = ["block_id", "clade", "family", "decoy", "reference_id",
              "n_references_contributing", "member_fraction",
              "n_members", "n_clade_windows",
              "recurrence", "coord_tightness", "mad_start_bp", "strand_consistency",
              "length_consistency", "composite_score", "median_start", "median_len",
              "median_pident", "ref_start", "ref_end", "ref_len", "dist_to_anchor",
              "side"]


STEP4_CFG = Path(__file__).resolve().parents[1] / "configs" / "blind" / "step4.json"


def load_cfg() -> dict:
    """Frozen defaults. The score floors are properties of the decoy null, so
    they travel with the other parameters or they mean nothing."""
    if STEP4_CFG.exists():
        return json.loads(STEP4_CFG.read_text())
    return {}


def main() -> int:
    cfg = load_cfg()
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", required=True, type=Path)
    ap.add_argument("--decoys", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--allow-partial", action="store_true",
                    help="keep output when some clades failed, for debugging; "
                         "the run is marked incomplete and candidate_table.py "
                         "will refuse it")
    ap.add_argument("--threads", type=int, default=4, help="blastn threads per clade")
    ap.add_argument("--references", type=int, default=cfg.get("k_references"))
    ap.add_argument("--proteins", type=Path,
                    help="dir of per-clade anchor-protein FASTAs; reference\n                          diversity is measured on protein, not on RNA")
    ap.add_argument("--evalue", type=float, default=cfg.get("evalue"))
    ap.add_argument("--min-len", type=int, default=cfg.get("min_len"))
    ap.add_argument("--min-recurrence", type=float, default=cfg.get("min_recurrence"))
    ap.add_argument("--min-members", type=int, default=cfg.get("min_members"))
    ap.add_argument("--max-coding-frac", type=float, default=cfg.get("max_coding_frac"))
    ap.add_argument("--tmp", type=Path, default=Path(os.environ.get("TMPDIR", "/tmp")))
    args = ap.parse_args()
    require_slurm("step 4 block discovery")
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "blocks").mkdir(exist_ok=True)

    if args.proteins and args.proteins.is_dir():
        for fp in sorted(args.proteins.glob("*.faa")):
            for h, s in read_fasta(str(fp)):
                _PROT.setdefault(h.split()[0], s)
        print(f"loaded {len(_PROT)} anchor proteins for reference diversity")
    else:
        print("no --proteins given; reference diversity falls back to window "
              "sequence, which is a weaker proxy", file=sys.stderr)

    tasks = [(p, False, args) for p in sorted(args.windows.glob("*.fasta"))]
    if args.decoys:
        tasks += [(p, True, args) for p in sorted(args.decoys.glob("*.fasta"))]
    print(f"{sum(1 for t in tasks if not t[1])} real clades, "
          f"{sum(1 for t in tasks if t[1])} decoy clades")

    all_blocks, all_stats, fails = [], [], []
    ctx = mp.get_context("fork")
    with ctx.Pool(args.workers) as pool:
        for i, (blocks, stats, err) in enumerate(
                pool.imap_unordered(_run, tasks), 1):
            if err:
                fails.append((stats["clade"], err))
            all_blocks.extend(blocks)
            all_stats.append(stats)
            if i % 20 == 0 or i == len(tasks):
                print(f"  {i}/{len(tasks)} clades, {len(all_blocks)} blocks",
                      flush=True)

    with open(args.out / "blocks.tsv", "w") as fh:
        fh.write("\t".join(BLOCK_COLS) + "\n")
        for b in sorted(all_blocks, key=lambda b: (b["decoy"], -b["composite_score"])):
            fh.write("\t".join(str(b[c]) for c in BLOCK_COLS) + "\n")

    # Per-anchor block coordinates -- what step 5-8 and the benchmark consume.
    with open(args.out / "block_members.tsv", "w") as fh:
        fh.write("block_id\tclade\tdecoy\tanchor_id\tstart\tend\tlen\tstrand\tpident\n")
        for b in all_blocks:
            for m in b["members"]:
                fh.write(f"{b['block_id']}\t{b['clade']}\t{b['decoy']}\t"
                         f"{m['anchor_id']}\t{m['start']}\t{m['end']}\t{m['len']}\t"
                         f"{m['strand']}\t{m['pident']}\n")

    # Per-block FASTA, real clades only -- decoy sequence has no downstream use.
    n_fa = 0
    cached_clade, recs = None, {}
    # Grouped by clade so each window FASTA is read once, not once per block.
    for b in sorted((b for b in all_blocks if not b["decoy"]),
                    key=lambda b: (b["clade"], -b["composite_score"])):
        if b["clade"] != cached_clade:
            cached_clade = b["clade"]
            recs = {r["anchor_id"]: r
                    for r in load_clade(args.windows / f"{b['clade']}.fasta")}
        with open(args.out / "blocks" / f"{b['block_id']}.fasta", "w") as fh:
            for m in b["members"]:
                r = recs.get(m["anchor_id"])
                if not r:
                    continue
                s = max(0, m["start"] + r["anchor_offset"])
                e = min(r["win_len"], m["end"] + r["anchor_offset"])
                if e - s < 10:
                    continue
                fh.write(f">{m['anchor_id']} block={b['block_id']} "
                         f"anchor_frame={m['start']}..{m['end']} "
                         f"strand={m['strand']}\n{r['seq'][s:e]}\n")
        n_fa += 1

    scols = ["clade", "decoy", "n_windows", "n_reference_windows", "references_used",
             "n_members_with_any_hsp", "frac_members_with_any_hsp",
             "n_members_with_noncoding_hsp", "frac_members_with_noncoding_hsp",
             "n_blocks", "top_block_id", "top_block_member_fraction"]
    with open(args.out / "clade_stats.tsv", "w") as fh:
        fh.write("\t".join(scols) + "\n")
        for s in sorted(all_stats, key=lambda s: (s.get("decoy"), s["clade"])):
            fh.write("\t".join(str(s.get(c, "")) for c in scols) + "\n")
    with open(args.out / "reference_recovery.tsv", "w") as fh:
        fh.write("clade\tdecoy\treference_id\trecovered_members\tn_windows\n")
        for s in all_stats:
            for rid, nrec in (s.get("per_reference_recovered_members") or {}).items():
                fh.write(f"{s['clade']}\t{s.get('decoy')}\t{rid}\t{nrec}\t"
                         f"{s.get('n_windows','')}\n")

    real = [b for b in all_blocks if not b["decoy"]]
    decoy = [b for b in all_blocks if b["decoy"]]
    summary = {
        "real_clades": sum(1 for s in all_stats if not s["decoy"]),
        "decoy_clades": sum(1 for s in all_stats if s["decoy"]),
        "real_blocks": len(real),
        "decoy_blocks": len(decoy),
        "real_clades_with_a_block": len({b["clade"] for b in real}),
        "decoy_clades_with_a_block": len({b["clade"] for b in decoy}),
        # The run is usable downstream only if every clade was processed. A
        # failed decoy clade silently shrinks the null's denominator and makes
        # the empirical FDR look better than it is; a failed real clade costs
        # sensitivity without saying which clade went missing. Neither shows up
        # in real_clades/decoy_clades, which count attempted clades, not
        # successful ones.
        "complete": not fails,
        "n_failures": len(fails),
        "failures": fails[:10],
        "median_frac_members_with_noncoding_hsp_real": (statistics.median(
            [s["frac_members_with_noncoding_hsp"] for s in all_stats
             if not s["decoy"] and "frac_members_with_noncoding_hsp" in s] or [0])),
        "median_n_reference_windows_real": (statistics.median(
            [s["n_reference_windows"] for s in all_stats
             if not s["decoy"] and "n_reference_windows" in s] or [0])),
        "params": {k: str(v) for k, v in vars(args).items() if k != "tmp"},
    }
    for label, bs in (("real", real), ("decoy", decoy)):
        if bs:
            sc = sorted(b["composite_score"] for b in bs)
            summary[f"{label}_composite"] = {
                "n": len(sc), "median": sc[len(sc) // 2],
                "p90": sc[int(len(sc) * 0.9)], "max": sc[-1]}
    (args.out / "discovery_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: v for k, v in summary.items() if k != "params"}, indent=2))
    print(f"{n_fa} block FASTAs -> {args.out}/blocks")
    print(f"-> {args.out}")

    if fails:
        with open(args.out / "failure_report.tsv", "w") as fh:
            fh.write("clade\tdecoy\terror\n")
            for clade, err in fails:
                fh.write(f"{clade}\t{clade.startswith('DECOY_')}\t{err}\n")
        n_d = sum(1 for c, _ in fails if c.startswith("DECOY_"))
        msg = (f"{len(fails)} of {len(tasks)} clades failed "
               f"({len(fails) - n_d} real, {n_d} decoy); "
               f"see {args.out / 'failure_report.tsv'}")
        if not args.allow_partial:
            # Exit non-zero so an sbatch chain stops here rather than building a
            # candidate table on an incomplete null.
            raise SystemExit(f"{msg}\nRefusing to report partial step-4 output. "
                             f"Re-run, or pass --allow-partial to keep it for "
                             f"debugging (it will be marked incomplete and "
                             f"candidate_table.py will refuse it).")
        print(f"WARNING: {msg}; output marked incomplete", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
