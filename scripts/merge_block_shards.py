#!/usr/bin/env python
"""Merge sharded discover_blocks outputs into one step-4 directory.

    merge_block_shards.py --shards <dir>/shard_* --out <dir>

Step 4 scores every clade independently: references are drawn from the clade's
own windows, and the anchor-protein pool handed to each shard is the same one.
Splitting the clade files across array tasks therefore produces exactly what one
job would have produced, and this puts the pieces back in the layout the floor,
null-report and lane-report steps already read.

What it refuses to do is paper over a partial run. If any shard reports a
failure or is missing its tables, that is carried into the merged summary and
the exit status, because a step-4 directory that is quietly short a few hundred
clades looks identical to a complete one and every downstream count -- the
floor, the recurrence denominators, the empirical bound -- would be computed on
it without complaint.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

TABLES = ["blocks.tsv", "block_members.tsv", "clade_stats.tsv",
          "reference_recovery.tsv"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shards", required=True, type=Path, nargs="+")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--expect-clades", type=int, default=0,
                    help="total real+decoy clades the shards should cover; a "
                         "mismatch is an error, not a warning")
    args = ap.parse_args()
    shards = [s for s in args.shards if s.is_dir()]
    if not shards:
        print("no shard directories found", file=sys.stderr)
        return 2
    (args.out / "blocks").mkdir(parents=True, exist_ok=True)

    problems = []
    for t in TABLES:
        header, n = None, 0
        with open(args.out / t, "w") as dst:
            for s in sorted(shards):
                p = s / t
                if not p.exists():
                    problems.append(f"{s.name}: missing {t}")
                    continue
                with open(p) as fh:
                    h = fh.readline()
                    if header is None:
                        header = h
                        dst.write(h)
                    elif h != header:
                        problems.append(f"{s.name}: {t} header differs")
                        continue
                    for line in fh:
                        if line.strip():
                            dst.write(line)
                            n += 1
        print(f"  {t}: {n} rows")

    nfa = 0
    for s in sorted(shards):
        d = s / "blocks"
        if d.is_dir():
            for f in d.glob("*.fasta"):
                shutil.copy(f, args.out / "blocks" / f.name)
                nfa += 1
    print(f"  blocks/: {nfa} block FASTAs")

    agg = {"real_clades": 0, "decoy_clades": 0, "real_blocks": 0,
           "decoy_blocks": 0, "real_clades_with_a_block": 0,
           "decoy_clades_with_a_block": 0, "n_failures": 0}
    failures, params, incomplete = [], None, []
    for s in sorted(shards):
        p = s / "discovery_summary.json"
        if not p.exists():
            problems.append(f"{s.name}: missing discovery_summary.json")
            continue
        d = json.loads(p.read_text())
        for k in agg:
            agg[k] += int(d.get(k, 0) or 0)
        failures += list(d.get("failures") or [])
        if not d.get("complete", True):
            incomplete.append(s.name)
        if params is None:
            params = d.get("params")
        elif params != d.get("params"):
            problems.append(f"{s.name}: step-4 params differ from shard 0")

    total = agg["real_clades"] + agg["decoy_clades"]
    if args.expect_clades and total != args.expect_clades:
        problems.append(f"expected {args.expect_clades} clades, merged {total}")

    agg.update({
        "complete": not problems and not incomplete and not failures,
        "failures": failures, "incomplete_shards": incomplete,
        "n_shards": len(shards), "params": params,
        "merge_problems": problems,
    })
    (args.out / "discovery_summary.json").write_text(json.dumps(agg, indent=2))
    print(f"\n{agg['real_clades']} real + {agg['decoy_clades']} decoy clades, "
          f"{agg['real_blocks']} real + {agg['decoy_blocks']} decoy blocks")
    if problems or failures or incomplete:
        for p in problems + [f"failure: {x}" for x in failures] + [
                f"incomplete shard: {x}" for x in incomplete]:
            print(f"  PROBLEM {p}", file=sys.stderr)
        return 1
    print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
