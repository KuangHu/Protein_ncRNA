#!/usr/bin/env python
"""Integrity report for genomes_db, and the canonical real-genome manifest.

74% of `INDEX.tsv` is 20-byte gzip stubs left by failed downloads, so the raw
row count overstates the mirror by ~6x. Anything that samples INDEX.tsv naively
gets a sample that is mostly dead files: effective n far below the nominal n,
species depth overstated, downstream hit rates divided by the wrong denominator,
and no clustering power because the same few real files keep being drawn.

`out/real_genomes.tsv` written here is meant to be the **shared entry point**
for every pipeline that reads this mirror. Prefer it over scanning INDEX.tsv.

    genomes_db_report.py --out out                   # size check only (fast)
    genomes_db_report.py --out out --verify-gzip     # also run gzip -t (slow)

`--verify-gzip` tests that each candidate decompresses, which catches truncated
downloads that a size check passes. It reads every file, so it is a SLURM job,
not a login-node one.
"""

from __future__ import annotations

import argparse
import gzip
import json
import multiprocessing as mp
import sys
from pathlib import Path

DB = Path("/global/scratch/users/kh36969/genomes_db")
MIN_BYTES = 100_000


def read_index(db: Path) -> list[tuple[str, str, int]]:
    rows = []
    with open(db / "INDEX.tsv") as fh:
        header = next(fh).rstrip("\n").split("\t")
        try:
            i_acc, i_path, i_bytes = (header.index(c) for c in ("accession", "rel_path", "bytes"))
        except ValueError:
            i_acc, i_path, i_bytes = 0, 1, 2
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_acc, i_path, i_bytes):
                continue
            n = int(f[i_bytes]) if f[i_bytes].isdigit() else -1
            rows.append((f[i_acc], f[i_path], n))
    return rows


def species_map(db: Path) -> dict[str, str]:
    """Accession -> species slug, from the per-species symlink views.

    The views hold ~611k symlinks, so this reads directory entries by name via
    scandir and never stats them; stat'ing each one takes minutes on this
    filesystem and tells us nothing we need.
    """
    import os

    out: dict[str, str] = {}
    views = db / "views"
    if not views.is_dir():
        return out
    for view in os.scandir(views):
        by_sp = Path(view.path) / "by_species"
        if not by_sp.is_dir():
            continue
        for sp in os.scandir(by_sp):
            if not sp.is_dir():
                continue
            for e in os.scandir(sp.path):
                if e.name.endswith(".fna.gz"):
                    out.setdefault(e.name[:-len(".fna.gz")], sp.name)
    return out


def _gzip_ok(task: tuple[str, str]) -> tuple[str, bool, int]:
    """Decompress a little of the file: catches truncated and corrupt downloads."""
    acc, path = task
    try:
        total = 0
        with gzip.open(path, "rb") as fh:
            while True:
                chunk = fh.read(1 << 20)
                if not chunk:
                    break
                total += len(chunk)
        return acc, total > 0, total
    except Exception:
        return acc, False, 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", type=Path, default=DB)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--min-bytes", type=int, default=MIN_BYTES)
    ap.add_argument("--check-exists", action="store_true",
                    help="stat every indexed path (slow on this filesystem)")
    ap.add_argument("--verify-gzip", action="store_true")
    ap.add_argument("--workers", type=int, default=16)
    args = ap.parse_args()
    if args.verify_gzip:
        require_slurm("genomes_db gzip verification")
    args.out.mkdir(parents=True, exist_ok=True)

    rows = read_index(args.db)
    sp_of = species_map(args.db)
    print(f"INDEX.tsv rows: {len(rows)}")

    # Stat'ing all ~305k paths costs minutes here and only confirms what
    # INDEX.tsv already records, so it is opt-in. --verify-gzip subsumes it for
    # the files that matter, since a missing file fails to open.
    present, big, missing_path = [], [], 0
    for acc, rel, n in rows:
        p = args.db / rel
        if args.check_exists and not p.exists():
            missing_path += 1
            continue
        present.append((acc, str(p), n))
        if n >= args.min_bytes:
            big.append((acc, str(p), n))
    if args.check_exists:
        print(f"files present on disk: {len(present)}   (missing: {missing_path})")
    else:
        print(f"index rows with a path: {len(present)}   (existence not checked; "
              f"use --check-exists)")
    print(f"files >= {args.min_bytes} B: {len(big)}")

    verified: dict[str, int] | None = None
    if args.verify_gzip:
        print(f"verifying gzip integrity of {len(big)} files on {args.workers} workers ...")
        verified = {}
        bad = 0
        ctx = mp.get_context("fork")
        with ctx.Pool(args.workers) as pool:
            for i, (acc, ok, nbytes) in enumerate(
                    pool.imap_unordered(_gzip_ok, [(a, p) for a, p, _ in big], chunksize=16), 1):
                if ok:
                    verified[acc] = nbytes
                else:
                    bad += 1
                if i % 5000 == 0:
                    print(f"  {i}/{len(big)} checked, {bad} bad", flush=True)
        print(f"gzip -t equivalent: {len(verified)} pass, {bad} fail")

    real = [(a, p, n) for a, p, n in big if verified is None or a in verified]

    per_species: dict[str, dict[str, int]] = {}
    for acc, _, _ in rows:
        sp = sp_of.get(acc, "(not in any view)")
        per_species.setdefault(sp, {"indexed": 0, "real": 0})["indexed"] += 1
    real_accs = {a for a, _, _ in real}
    for acc in real_accs:
        sp = sp_of.get(acc, "(not in any view)")
        per_species.setdefault(sp, {"indexed": 0, "real": 0})["real"] += 1

    manifest = args.out / "real_genomes.tsv"
    with open(manifest, "w") as fh:
        fh.write("accession\tpath\tbytes\tspecies\n")
        for acc, p, n in sorted(real):
            fh.write(f"{acc}\t{p}\t{n}\t{sp_of.get(acc, 'unknown')}\n")

    print(f"\n{'species':<32}{'indexed':>10}{'real':>10}{'%real':>8}")
    for sp, d in sorted(per_species.items(), key=lambda kv: -kv[1]["real"]):
        pct = 100 * d["real"] / d["indexed"] if d["indexed"] else 0.0
        print(f"{sp:<32}{d['indexed']:>10}{d['real']:>10}{pct:>7.1f}%")

    summary = {
        "index_rows": len(rows),
        "files_present": len(present) if args.check_exists else None,
        "files_missing": missing_path if args.check_exists else None,
        "files_ge_min_bytes": len(big),
        "min_bytes": args.min_bytes,
        "gzip_verified": len(verified) if verified is not None else None,
        "real_genomes": len(real),
        "pct_real": round(100 * len(real) / max(len(rows), 1), 2),
        "per_species": per_species,
        "manifest": str(manifest),
    }
    (args.out / "genomes_db_report.json").write_text(json.dumps(summary, indent=2))
    print(f"\nreal genomes: {len(real)} / {len(rows)} index rows "
          f"({summary['pct_real']}%)\nmanifest -> {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
