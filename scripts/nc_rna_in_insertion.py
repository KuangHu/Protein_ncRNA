#!/usr/bin/env python
"""Find ncRNA in the non-coding regions of the insertion. Nothing else.

    nc_rna_in_insertion.py --sites <insertions_sites.jsonl> --out <dir> \\
                           [--bags <list>] [--min-recurrence 0.4] \\
                           [--min-sites 5] [--decoys <n>]

The insertion is `[nc][tnp][nc]`. The insertion caller already resolved it: it
compared the occupied locus to the empty one to define the element, called the
CDS inside, and wrote what is left as `noncoding_regions` -- per site, already
in-element, already non-coding, median 116 bp.

So the search space is those regions and only those regions. There is no host
flank in it, no side to choose and no element extent to estimate.

A candidate is a non-coding region that RECURS: the same sequence appearing in
independent insertion events of the same protein family, in different genomes.
One site shows nothing; a region conserved across many independent insertions of
the same transposase is under selection, and that is the signal.

Recurrence is counted over DISTINCT SITES, never over sequences. One site
contributing three near-identical regions is one observation, and counting it as
three would manufacture recurrence out of a single insertion.

The null permutes regions between bags of the same family, keeping cluster sizes
and the bag's site count, so a decoy faces the same question with the same
amount of data. Nothing in here reads a known guide RNA, msr/msd, tracrRNA or
bridge RNA sequence, length or structure.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import subprocess
import sys
from collections import defaultdict
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


def cluster(seqs, tmp: Path, mmseqs: str, ident: float, cov: float,
            threads: int):
    """Cluster (site_id, seq) pairs; return list of clusters of site_ids.

    Identical sequences are collapsed first so a bag with thousands of copies of
    one region does not pay for them in the all-vs-all.
    """
    uniq = {}
    for sid, s in seqs:
        uniq.setdefault(s, []).append(sid)
    if len(uniq) < 2:
        return [list({x for v in uniq.values() for x in v})] if uniq else []
    fa = tmp / "in.fa"
    with open(fa, "w") as fh:
        for i, s in enumerate(uniq):
            fh.write(f">s{i}\n{s}\n")
    keys = list(uniq)
    out = tmp / "clu"
    r = subprocess.run(
        [mmseqs, "easy-cluster", str(fa), str(out), str(tmp / "mm"),
         "--min-seq-id", str(ident), "-c", str(cov), "--cov-mode", "1",
         "--threads", str(threads), "-v", "0"],
        capture_output=True, text=True)
    tsv = Path(str(out) + "_cluster.tsv")
    if r.returncode != 0 or not tsv.exists():
        return []
    groups = defaultdict(set)
    for line in open(tsv):
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        rep, mem = f[0], f[1]
        for sid in uniq[keys[int(mem[1:])]]:
            groups[rep].add(sid)
    return [sorted(v) for v in groups.values()]


def best_for(seqs, n_sites, tmp, mmseqs, args, want_members=False):
    """Strongest cluster of a clade: (recurrence, n_sites, median_len[, seqs]).

    With want_members the cluster's own sequences come back too. Stage 2 counts
    this candidate in the tnp-anchored windows, which are 7x more numerous than
    the called insertion sites, so the sequence has to leave this script or that
    stage has nothing to search for.
    """
    cl = cluster(seqs, tmp, mmseqs, args.ident, args.cov, args.threads)
    by_site = defaultdict(list)
    for sid, s in seqs:
        by_site[sid].append(s)
    best = None
    for g in cl:
        k = len(set(g))
        if k < args.min_sites:
            continue
        rec = k / n_sites
        mem = [s for sid in g for s in by_site.get(sid, [])]
        lens = sorted(len(s) for s in mem)
        if best is None or rec > best[0]:
            med = lens[len(lens) // 2] if lens else 0
            best = (rec, k, med, mem) if want_members else (rec, k, med)
    return best


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--bags", type=Path)
    ap.add_argument("--min-len", type=int, default=50)
    ap.add_argument("--min-sites", type=int, default=5,
                    help="distinct insertion events a candidate must appear in")
    ap.add_argument("--min-recurrence", type=float, default=0.4)
    ap.add_argument("--ident", type=float, default=0.70)
    ap.add_argument("--cov", type=float, default=0.70)
    ap.add_argument("--decoys", type=int, default=3,
                    help="same-family permutation decoys per bag")
    ap.add_argument("--max-sites", type=int, default=3000)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1704)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    args = ap.parse_args()
    require_slurm()
    args.out.mkdir(parents=True, exist_ok=True)
    mmseqs = "mmseqs"

    want = set(args.bags.read_text().split()) if args.bags else None
    # bag -> site_id -> [nc region sequences]
    bags = defaultdict(lambda: defaultdict(list))
    for line in open(args.sites):
        try:
            r = json.loads(line)
        except json.JSONDecodeError:
            continue
        b = r.get("bag_id")
        if not b or (want is not None and b not in want):
            continue
        for s in (r.get("noncoding_regions") or []):
            if len(s) >= args.min_len:
                bags[b][r["site_id"]].append(s)
    allb = sorted(b for b in bags if len(bags[b]) >= args.min_sites)
    log(f"{len(allb)} bags with >= {args.min_sites} sites carrying "
        f"non-coding >= {args.min_len} bp")

    mine = [b for i, b in enumerate(allb) if i % args.nshards == args.shard]
    rng = random.Random(args.seed)
    tmp = args.out / f"_tmp{args.shard}"
    tmp.mkdir(exist_ok=True)
    cdir = args.out / "candidates"
    cdir.mkdir(exist_ok=True)

    # family -> bags, for the permutation null
    byfam = defaultdict(list)
    for b in allb:
        byfam[family_of(b)].append(b)

    rows = []
    for bi, bag in enumerate(mine, 1):
        sites = bags[bag]
        sids = sorted(sites)
        if len(sids) > args.max_sites:
            sids = rng.sample(sids, args.max_sites)
        seqs = [(sid, s) for sid in sids for s in sites[sid]]
        n = len(sids)
        real = best_for(seqs, n, tmp, mmseqs, args, want_members=True)

        # same-family permutation: draw regions from OTHER bags of the family,
        # matched on this bag's site count, so the decoy answers the same
        # question with the same amount of data.
        pool = [x for ob in byfam[family_of(bag)] if ob != bag
                for x in bags[ob].items()]
        dec = []
        for d in range(args.decoys if len(pool) >= n else 0):
            pick = rng.sample(pool, n)
            dseqs = [(f"d{d}_{i}", s) for i, (sid, ss) in enumerate(pick)
                     for s in ss]
            db = best_for(dseqs, n, tmp, mmseqs, args)
            if db:
                dec.append(db[0])
        # The candidate's own sequences, for stage 2. A median-length member
        # is the representative: a cluster consensus would be a sequence that
        # exists in no genome, and stage 2 searches real windows.
        if real:
            mem = sorted(real[3], key=len)
            cand = mem[len(mem) // 2]
            with open(cdir / f"{bag}.fasta", "w") as fh:
                fh.write(f">{bag}|candidate|rec={real[0]:.4f}|"
                         f"sites={real[1]}|len={len(cand)}\n{cand}\n")
                for i, m in enumerate(mem[:200]):
                    fh.write(f">{bag}|member{i}|len={len(m)}\n{m}\n")
        rows.append({
            "bag": bag, "n_sites": n,
            "n_nc_regions": len(seqs),
            "recurrence": round(real[0], 4) if real else 0.0,
            "n_sites_with_candidate": real[1] if real else 0,
            "candidate_len": real[2] if real else 0,
            "n_decoys": len(dec),
            "decoy_max_recurrence": round(max(dec), 4) if dec else "",
            "margin": round(real[0] - max(dec), 4) if (real and dec) else "",
            "passes": bool(real and real[0] >= args.min_recurrence
                           and (not dec or real[0] > max(dec))),
        })
        if bi % 25 == 0 or bi == len(mine):
            log(f"  [{bi}/{len(mine)}] {bag}: "
                f"rec {rows[-1]['recurrence']} "
                f"decoy {rows[-1]['decoy_max_recurrence']}")

    cols = list(rows[0]) if rows else []
    p = args.out / f"nc_candidates_{args.shard:02d}.tsv"
    with open(p, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    npass = sum(1 for r in rows if r["passes"])
    log(f"\n{len(rows)} bags, {npass} with a recurrent in-element "
        f"non-coding candidate that beat its null")
    log(f"-> {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
