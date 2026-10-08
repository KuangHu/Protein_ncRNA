#!/usr/bin/env python
"""Assemble an Arm 2 census in the layout extract_windows.py expects.

    build_census.py --anchors <expand/anchors> --bags <bag list> --out <dir>

The sharded anchor table and per-shard protein FASTAs are collapsed into the
one `anchors.tsv` + `proteins/<bag>.faa` layout the window extractor reads. The
bag list must already include the sibling bags of every family the permutation
null will use: a sibling with no windows is silently skipped downstream, which
turns into a missing null rather than an error.
"""

from __future__ import annotations

import argparse
import glob
import shutil
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--anchors", required=True, type=Path)
    ap.add_argument("--bags", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    (args.out / "proteins").mkdir(parents=True, exist_ok=True)
    want = set(args.bags.read_text().split())

    hdr, n = None, 0
    with open(args.out / "anchors.tsv", "w") as out:
        for p in sorted(glob.glob(str(args.anchors / "anchors_shard*.tsv"))):
            with open(p) as fh:
                h = next(fh)
                if hdr is None:
                    hdr = h
                    out.write(h)
                    ci = {k: i for i, k in
                          enumerate(h.rstrip("\n").split("\t"))}
                # Split only as far as the family column; these shards are wide
                # and the full split is most of the runtime.
                for line in fh:
                    if line.split("\t", ci["family"] + 1)[ci["family"]] in want:
                        out.write(line)
                        n += 1
    print(f"{n} anchor rows for {len(want)} bags", flush=True)

    missing = []
    for bag in sorted(want):
        parts = sorted(glob.glob(
            str(args.anchors / "proteins" / f"{bag}.shard*.faa")))
        if not parts:
            missing.append(bag)
            continue
        with open(args.out / "proteins" / f"{bag}.faa", "wb") as out:
            for p in parts:
                with open(p, "rb") as fh:
                    shutil.copyfileobj(fh, out)
    print(f"{len(want) - len(missing)} protein FASTAs written")
    if missing:
        # Loud, because a bag with no proteins falls back to DNA references and
        # would then be the only arm in a different reference space.
        print(f"WARNING {len(missing)} bags have no protein shards: "
              f"{missing[:10]}")
    print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
