#!/usr/bin/env python
"""Trim each bag's windows to the insertion, before anything is scored.

    trim_windows_to_element.py --geometry <element_geometry.tsv> \\
                               --windows <dir ...> --out <dir> [--min-len 300]

The insertion filter was applied after scoring, which showed that 388 of 579
blocks were host flank -- but the floor, the decoy distribution and the whole
ranking had already been computed over that population. A post-hoc filter cannot
fix a null drawn from the wrong space.

This moves the constraint upstream. Each window is cut to
`[anchor - upstream_extent, anchor_end + downstream_extent]`, with each side
opened only where element_geometry found the element extends, so a `[tnp][nc]`
bag keeps nothing upstream of its tnp. Decoys built from these windows are drawn
from the same space the real clades occupy, which is the only way the floor
means anything.

anchor_offset is rewritten to the new frame. Getting that wrong would move every
block relative to its anchor while leaving every length intact -- the same
silent class of error as the plus-strand window bug -- so the anchor's own
sequence is re-checked after the cut and a window that fails is dropped rather
than emitted.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry", required=True, type=Path)
    ap.add_argument("--windows", required=True, type=Path, nargs="+")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--min-len", type=int, default=300,
                    help="a trimmed window shorter than this cannot hold a "
                         "block plus its anchor and is dropped")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    args = ap.parse_args()

    geo = {r["bag"]: r for r in csv.DictReader(
        open(args.geometry), delimiter="\t")}
    src = {}
    for root in args.windows:
        if root.is_dir():
            for d in sorted(root.iterdir()):
                if d.is_dir() and (d / f"{d.name}.jsonl").exists():
                    src.setdefault(d.name, d)

    bags = [b for i, b in enumerate(sorted(geo)) if i % args.nshards ==
            args.shard]
    args.out.mkdir(parents=True, exist_ok=True)
    tot_in = tot_out = dropped = 0
    for bag in bags:
        d = src.get(bag)
        if d is None:
            continue
        g = geo[bag]
        up = int(g["upstream_extent_bp"])
        dn = int(g["downstream_extent_bp"])
        od = args.out / bag
        (od / "proteins").mkdir(parents=True, exist_ok=True)
        kept = []
        with open(d / f"{bag}.jsonl") as fh:
            for line in fh:
                w = json.loads(line)
                tot_in += 1
                ao, al, seq = w["anchor_offset"], w["anchor_len"], w["seq"]
                s0 = max(0, ao - up)
                s1 = min(len(seq), ao + al + dn)
                if s1 - s0 < args.min_len:
                    dropped += 1
                    continue
                sub = seq[s0:s1]
                nao = ao - s0
                # The anchor must still be the anchor after the cut.
                if sub[nao:nao + al] != seq[ao:ao + al]:
                    dropped += 1
                    continue
                w2 = dict(w)
                w2["seq"] = sub
                w2["anchor_offset"] = nao
                w2["trimmed_from"] = [s0, s1]
                w2["element_upstream_bp"] = up
                w2["element_downstream_bp"] = dn
                kept.append(w2)
        if not kept:
            continue
        with open(od / f"{bag}.jsonl", "w") as fh:
            for w in kept:
                fh.write(json.dumps(w) + "\n")
        with open(od / f"{bag}.fasta", "w") as fh:
            for w in kept:
                fh.write(f">{w['anchor_id']} family={bag} "
                         f"anchor_offset={w['anchor_offset']}\n{w['seq']}\n")
        p = d / "proteins" / f"{bag}.faa"
        if p.exists():
            (od / "proteins" / f"{bag}.faa").write_bytes(p.read_bytes())
        tot_out += len(kept)
        print(f"  {bag}: {len(kept)} windows, up {up} dn {dn}, "
              f"median len {sorted(len(w['seq']) for w in kept)[len(kept)//2]}")
    print(f"\n{len(bags)} bags: {tot_in} windows in, {tot_out} out, "
          f"{dropped} dropped")
    print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
