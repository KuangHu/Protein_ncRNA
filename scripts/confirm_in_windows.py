#!/usr/bin/env python
"""Stage 2: count the insertion-derived candidate in the tnp-anchored windows.

    confirm_in_windows.py --candidates <dir> --windows <dir ...> \\
                          --out <dir> [--min-identity 80] [--min-cov 0.5]

Stage 1 defines a candidate inside the insertion, where the element boundary is
known, and it pays for that with sample size: a bag has a median of 17 called
insertion sites. Blasting the transposase back into the genomes finds a median
of 122 copies of it -- 7.1x more -- but a tnp hit gives no element boundary, and
carving a window around it is what produced two thirds host flank before.

The boundary is only needed to FIND a candidate, not to COUNT one. Asking
whether a known 191 bp sequence is present beside this tnp copy, on the same
side and in the same orientation, is a sequence question, and it can use every
window without inheriting the flank problem. So: boundary from the insertions,
statistics from the windows.

What that buys, beyond a larger denominator:

  FIXED INSERTIONS. The insertion caller only sees polymorphic sites -- ones
  with an observed empty allele somewhere in the sample. A transposase at a
  fixed insertion is invisible to it. A candidate that also appears at those
  copies is a component of the element, not an artifact of polymorphic sites.

  AN OUT-OF-SAMPLE TEST. The candidate was defined on one set of loci and is
  counted on a larger, mostly disjoint one.

Windows are already oriented: extract_windows.py normalises `anchor_strand`, so
a plus-strand hit is the same direction as the transposase. Side and orientation
are therefore read, not assumed, and a candidate whose hits scatter across both
sides and both strands is reported as such -- consistency is evidence, and its
absence has to be visible rather than averaged away.

The null is the same candidate searched in the windows of OTHER bags of the same
family, where the element is not. Nothing here re-enters stage 1.
"""

from __future__ import annotations

import argparse
import json
import os
import statistics as st
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


def log(*a):
    print(*a, flush=True)


def require_slurm():
    if os.environ.get("PN_ALLOW_LOGIN_NODE") == "1":
        return
    if not os.environ.get("SLURM_JOB_ID"):
        sys.exit("refusing to run on a login node; submit with sbatch")


def family_of(name: str) -> str:
    return name.split("__")[0]


def read_candidate(p: Path):
    """The representative only -- the first record."""
    name, buf = None, []
    for line in open(p):
        if line.startswith(">"):
            if name:
                break
            name = line[1:].strip()
        else:
            buf.append(line.strip())
    return "".join(buf)


def window_fasta(d: Path, bag: str, tmp: Path, cap: int):
    """Write a bag's windows, keeping each window's anchor offset and length."""
    meta = {}
    out = tmp / f"{bag}.win.fa"
    n = 0
    with open(out, "w") as fh:
        with open(d / f"{bag}.jsonl") as src:
            for line in src:
                if n >= cap:
                    break
                w = json.loads(line)
                wid = f"w{n}"
                meta[wid] = (int(w["anchor_offset"]), int(w["anchor_len"]))
                fh.write(f">{wid}\n{w['seq']}\n")
                n += 1
    return out, meta


def search(cand: str, subj: Path, tmp: Path, threads: int):
    q = tmp / "cand.fa"
    q.write_text(f">cand\n{cand}\n")
    r = subprocess.run(
        ["blastn", "-task", "blastn", "-word_size", "11", "-dust", "no",
         "-evalue", "1e-5", "-num_threads", str(threads),
         "-max_target_seqs", "100000", "-query", str(q), "-subject", str(subj),
         "-outfmt", "6 sseqid pident length sstart send"],
        capture_output=True, text=True)
    return r.stdout


def best_per_window(raw: str, clen: int, min_id: float, min_cov: float):
    best = {}
    for line in raw.splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        sid, pid, alen = f[0], float(f[1]), int(f[2])
        ss, se = int(f[3]), int(f[4])
        if pid < min_id or alen < min_cov * clen:
            continue
        if sid not in best or alen > best[sid][1]:
            best[sid] = (pid, alen, ss, se)
    return best


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True, type=Path)
    ap.add_argument("--windows", required=True, type=Path, nargs="+")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--min-identity", type=float, default=80.0)
    ap.add_argument("--min-cov", type=float, default=0.5)
    ap.add_argument("--max-windows", type=int, default=4000)
    ap.add_argument("--n-offbag", type=int, default=3,
                    help="same-family bags searched as the null")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    args = ap.parse_args()
    require_slurm()
    args.out.mkdir(parents=True, exist_ok=True)

    src = {}
    for root in args.windows:
        if root.is_dir():
            for d in sorted(root.iterdir()):
                if d.is_dir() and (d / f"{d.name}.jsonl").exists():
                    src.setdefault(d.name, d)
    cands = {p.stem: p for p in sorted(args.candidates.glob("*.fasta"))}
    both = sorted(set(cands) & set(src))
    log(f"{len(cands)} stage-1 candidates, {len(src)} bags with windows, "
        f"{len(both)} testable")

    byfam = defaultdict(list)
    for b in src:
        byfam[family_of(b)].append(b)

    mine = [b for i, b in enumerate(both) if i % args.nshards == args.shard]
    tmp = Path(tempfile.mkdtemp())
    rows = []
    for bi, bag in enumerate(mine, 1):
        cand = read_candidate(cands[bag])
        if len(cand) < 40:
            continue
        subj, meta = window_fasta(src[bag], bag, tmp, args.max_windows)
        hits = best_per_window(search(cand, subj, tmp, args.threads),
                               len(cand), args.min_identity, args.min_cov)
        sides, strands, dists = Counter(), Counter(), []
        for wid, (pid, alen, ss, se) in hits.items():
            ao, al = meta[wid]
            plus = se > ss
            lo, hi = (ss, se) if plus else (se, ss)
            # anchor frame: the tnp occupies [0, anchor_len)
            a0, a1 = lo - 1 - ao, hi - ao
            side = "upstream" if a1 <= 0 else (
                "downstream" if a0 >= al else "overlaps_tnp")
            sides[side] += 1
            strands["same_as_tnp" if plus else "opposite"] += 1
            dists.append(abs(a0) if side == "upstream" else abs(a0 - al))

        # null: the same candidate in other bags of the same family
        fam = [b for b in byfam[family_of(bag)] if b != bag][:args.n_offbag]
        off_w = off_h = 0
        for ob in fam:
            osubj, ometa = window_fasta(src[ob], ob, tmp, args.max_windows)
            oh = best_per_window(search(cand, osubj, tmp, args.threads),
                                 len(cand), args.min_identity, args.min_cov)
            off_w += len(ometa)
            off_h += len(oh)
            osubj.unlink(missing_ok=True)

        nw = len(meta)
        modal = sides.most_common(1)[0] if sides else ("none", 0)
        rate = len(hits) / nw if nw else 0.0
        orate = off_h / off_w if off_w else 0.0
        rows.append({
            "bag": bag, "candidate_len": len(cand),
            "n_windows": nw, "n_windows_with_candidate": len(hits),
            "window_recurrence": round(rate, 4),
            "modal_side": modal[0],
            "frac_modal_side": round(modal[1] / len(hits), 3) if hits else 0.0,
            "frac_same_strand_as_tnp": round(
                strands["same_as_tnp"] / len(hits), 3) if hits else 0.0,
            "median_dist_to_tnp": int(st.median(dists)) if dists else "",
            "n_offbag_windows": off_w, "n_offbag_hits": off_h,
            "offbag_rate": round(orate, 4),
            "enrichment_over_offbag": round(rate / orate, 2) if orate else "",
        })
        subj.unlink(missing_ok=True)
        if bi % 20 == 0 or bi == len(mine):
            log(f"  [{bi}/{len(mine)}] {bag}: {len(hits)}/{nw} windows "
                f"({rate:.3f}), offbag {orate:.4f}")

    cols = list(rows[0]) if rows else []
    p = args.out / f"window_confirmation_{args.shard:02d}.tsv"
    with open(p, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    log(f"\n{len(rows)} candidates counted in tnp-anchored windows -> {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
