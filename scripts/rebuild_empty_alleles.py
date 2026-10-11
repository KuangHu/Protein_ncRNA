#!/usr/bin/env python
"""Re-extract each insertion's EMPTY allele wide, from the genomes themselves.

    rebuild_empty_alleles.py --sites <insertions_sites.jsonl> --bags <list> \\
                             --genome-map <tsv> --out <dir> [--flank 5000]

bags_v4 already pairs every occupied site with a real observed empty allele --
78,120 of them, none inferred. But `flank_len` is 120 and
`insertion_point_in_flank` is 60 for every single record, so what is stored is
60 bp of host on each side of the junction and nothing else.

That is too narrow to test a candidate. In-insertion blocks sit a median 458 bp
from the anchor (p90 1581), so not one of them reaches into the stored empty
flank. The occupied/empty pairs cannot currently confirm that any candidate
ncRNA is element-borne.

Widening them turns the pair into the decisive test, because the empty allele is
the same locus in a genome that never received the element:

    block present in occupied AND absent from empty  -> element-borne
    block present in BOTH                            -> host sequence

which is an observation, not an inference from window geometry. It settles the
79 `novel_in_element` blocks the annotation missed, and it independently retests
the 45 blocks rejected as beyond the element extent.

HOW AN EMPTY ALLELE IS RECOGNISED. The 120 bp flank is aligned to the species'
genomes and a hit counts only if ONE ungapped HSP spans the whole query. That is
the test itself, not a quality filter: a genome carrying the element cannot align
across the junction in one piece, because the element sits between the two
halves. So a full-length ungapped hit IS the absence of the insertion.

Nothing here reads a block, a score or a rank, so it cannot steer discovery; it
only rebuilds the locus the insertion caller already called.
"""

from __future__ import annotations

import argparse
import csv
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
        sys.exit("refusing to run on a login node; submit with sbatch "
                 "(or set PN_ALLOW_LOGIN_NODE=1 for a small deliberate test)")


def build_species_db(acc_paths, dest: Path, threads: int):
    """Concatenate genomes into one fasta, then makeblastdb + faidx on it.

    One file per species keeps the blast database and the extraction index in
    the same coordinate space, so a hit can be sliced without re-opening the
    genome it came from.
    """
    fa = dest.with_suffix(".fna")
    if not (fa.exists() and Path(str(fa) + ".fai").exists()
            and Path(str(dest) + ".nin").exists()):
        n = 0
        with open(fa, "w") as out:
            for acc, p in acc_paths:
                try:
                    with open(p) as fh:
                        for line in fh:
                            if line.startswith(">"):
                                sid = line[1:].split()[0]
                                out.write(f">{acc}|{sid}\n")
                            else:
                                out.write(line)
                    n += 1
                except OSError:
                    continue
        log(f"    concatenated {n} genomes -> {fa} "
            f"({fa.stat().st_size / 1e9:.1f} GB)")
        subprocess.run(["samtools", "faidx", str(fa)], check=True)
        subprocess.run(["makeblastdb", "-in", str(fa), "-dbtype", "nucl",
                        "-out", str(dest)], check=True,
                       stdout=subprocess.DEVNULL)
    return fa


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", required=True, type=Path)
    ap.add_argument("--bags", required=True, type=Path)
    ap.add_argument("--genome-map", required=True, type=Path,
                    help="tsv: species <tab> directory of .fna genomes")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--flank", type=int, default=5000,
                    help="bp of empty-allele sequence to keep each side of the "
                         "insertion point")
    ap.add_argument("--max-genomes", type=int, default=2000,
                    help="genomes sampled per species for the search database; "
                         "the median site has 9 empty carriers, so a sample "
                         "finds them without indexing 20k genomes")
    ap.add_argument("--max-carriers", type=int, default=5,
                    help="empty alleles kept per site")
    ap.add_argument("--min-identity", type=float, default=95.0)
    ap.add_argument("--seed", type=int, default=1704)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--species", nargs="*",
                    help="restrict to these species")
    args = ap.parse_args()
    require_slurm()
    args.out.mkdir(parents=True, exist_ok=True)
    dbdir = args.out / "db"
    dbdir.mkdir(exist_ok=True)

    gmap = {}
    for r in csv.reader(open(args.genome_map), delimiter="\t"):
        if len(r) >= 2 and not r[0].startswith("#"):
            gmap[r[0]] = Path(r[1])

    bags = set(args.bags.read_text().split())
    sites = defaultdict(list)
    for line in open(args.sites):
        try:
            r = json.loads(line)
        except json.JSONDecodeError:
            continue
        if r.get("bag_id") not in bags:
            continue
        if int(r.get("flank_len") or 0) != len(r.get("flank") or ""):
            continue
        sites[r.get("species")].append(r)
    log(f"{sum(len(v) for v in sites.values())} sites over {len(bags)} bags, "
        f"{len(sites)} species")

    rng = random.Random(args.seed)
    rows, nfa = [], 0
    fadir = args.out / "empty_alleles"
    fadir.mkdir(exist_ok=True)

    for sp in sorted(sites):
        if args.species and sp not in args.species:
            continue
        srcs = gmap.get(sp)
        if srcs is None or not srcs.is_dir():
            log(f"  {sp}: no genome directory in the map, skipped "
                f"({len(sites[sp])} sites)")
            for r in sites[sp]:
                rows.append({"site_id": r["site_id"], "bag": r["bag_id"],
                             "species": sp, "n_empty_alleles": 0,
                             "status": "no_genome_source"})
            continue
        gens = sorted(p for p in srcs.iterdir() if p.suffix == ".fna")
        if len(gens) > args.max_genomes:
            gens = rng.sample(gens, args.max_genomes)
        log(f"  {sp}: {len(sites[sp])} sites, {len(gens)} genomes sampled")
        fa = build_species_db([(p.stem, p) for p in sorted(gens)],
                              dbdir / sp, args.threads)

        q = args.out / f"{sp}.flanks.fa"
        with open(q, "w") as fh:
            for r in sites[sp]:
                fh.write(f">{r['site_id']}\n{r['flank']}\n")
        hits = subprocess.run(
            ["blastn", "-task", "blastn", "-word_size", "11", "-dust", "no",
             "-evalue", "1e-20", "-num_threads", str(args.threads),
             "-max_target_seqs", "200", "-db", str(dbdir / sp),
             "-query", str(q), "-outfmt",
             "6 qseqid sseqid pident length gapopen qstart qend sstart send"],
            capture_output=True, text=True).stdout

        per = defaultdict(list)
        for line in hits.splitlines():
            f = line.split("\t")
            if len(f) < 9:
                continue
            qid, sid = f[0], f[1]
            pid, alen, gaps = float(f[2]), int(f[3]), int(f[4])
            qs, qe, ss, se = (int(x) for x in f[5:9])
            # One ungapped HSP across the whole junction, or the element is
            # there and this is not an empty allele.
            if gaps or pid < args.min_identity or qs > 3 or qe < 118:
                continue
            per[qid].append((sid, qs, ss, se))

        idx = {r["site_id"]: r for r in sites[sp]}
        for sid_q, hs in per.items():
            rec = idx[sid_q]
            ip = int(rec.get("insertion_point_in_flank") or 60)
            out = []
            for sid, qs, ss, se in hs[:args.max_carriers]:
                plus = se > ss
                off = ip - qs            # bases from the HSP start to the point
                gp = ss + off if plus else ss - off
                lo, hi = gp - args.flank, gp + args.flank
                reg = f"{sid}:{max(1, lo)}-{hi}"
                seq = subprocess.run(["samtools", "faidx", str(fa), reg],
                                     capture_output=True, text=True).stdout
                body = "".join(seq.splitlines()[1:])
                if len(body) < 2 * args.flank * 0.25:
                    continue           # contig ended; too little to test against
                out.append((sid, reg, "+" if plus else "-", body))
            if out:
                with open(fadir / f"{sid_q}.fasta", "w") as fh:
                    for sid, reg, st, body in out:
                        fh.write(f">{sid_q}|{sid}|{reg}|{st}\n{body}\n")
                nfa += 1
            rows.append({"site_id": sid_q, "bag": rec["bag_id"], "species": sp,
                         "n_empty_alleles": len(out),
                         "status": "widened" if out else "no_empty_allele_found"})
        for r in sites[sp]:
            if r["site_id"] not in per:
                rows.append({"site_id": r["site_id"], "bag": r["bag_id"],
                             "species": sp, "n_empty_alleles": 0,
                             "status": "no_full_length_hit"})

    cols = ["site_id", "bag", "species", "n_empty_alleles", "status"]
    with open(args.out / "empty_allele_index.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    st = defaultdict(int)
    for r in rows:
        st[r["status"]] += 1
    log(f"\n{len(rows)} sites processed, {nfa} with a widened empty allele")
    for k, v in sorted(st.items(), key=lambda kv: -kv[1]):
        log(f"  {k:24s} {v:6d}  {100 * v / len(rows):5.1f}%")
    log(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
