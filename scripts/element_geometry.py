#!/usr/bin/env python
"""Per-bag element geometry: which side the non-coding sits on, and how far.

    element_geometry.py --bags <list> --windows <dir ...> \\
                        --sites <insertions_sites.jsonl> --out <dir> \\
                        [--sample 20] [--slack 50]

The insertion is `[non-coding][tnp]` or `[tnp][non-coding]`, and only non-coding
inside it counts. To restrict the search to the element we need two numbers per
bag: which side of the tnp the element extends to, and how far.

They are obtained differently on purpose.

  SIDE comes from sequence. The bag's recorded in-insertion non-coding regions
  are aligned to a sample of its windows and their offsets relative to the
  anchor are collected. Upstream hits mean `[non-coding][tnp]`, downstream mean
  `[tnp][non-coding]`, and a bag can legitimately show both.

  EXTENT does NOT come from sequence. It is `inserted_len - anchor_len`: the
  element's own measured length, from the insertion caller, minus the tnp. If
  the extent were taken from where the recorded non-coding aligned, the trimmed
  window would BE the recorded non-coding region and step 4 could only re-find
  what bags_v4 already annotated. Sizing by element length instead leaves the
  search blind inside the element, so a non-coding region the insertion caller
  missed can still be found.

So sequence decides where to look, and the element's length decides how much to
look at. Nothing here tells the block caller what to find.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics as st
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bags", required=True, type=Path)
    ap.add_argument("--windows", required=True, type=Path, nargs="+")
    ap.add_argument("--sites", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--sample", type=int, default=20,
                    help="windows per bag used to locate the recorded "
                         "non-coding; side is a property of the element, so a "
                         "sample settles it")
    ap.add_argument("--slack", type=int, default=50,
                    help="junction uncertainty added to each accepted side")
    ap.add_argument("--min-side-frac", type=float, default=0.15,
                    help="a side is accepted when it carries at least this "
                         "share of located non-coding, so a handful of stray "
                         "hits does not reopen the host flank")
    ap.add_argument("--max-nc", type=int, default=400)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    bags = set(args.bags.read_text().split())
    nc = defaultdict(dict)
    ins_len = defaultdict(list)
    with open(args.sites) as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            b = r.get("bag_id")
            if b not in bags:
                continue
            if r.get("inserted_len"):
                ins_len[b].append(int(r["inserted_len"]))
            if len(nc[b]) < args.max_nc:
                for s in (r.get("noncoding_regions") or []):
                    if len(s) >= 30:
                        nc[b].setdefault(s, None)

    src = {}
    for root in args.windows:
        if root.is_dir():
            for d in sorted(root.iterdir()):
                if d.is_dir() and (d / f"{d.name}.jsonl").exists():
                    src.setdefault(d.name, d)

    tmp = Path(tempfile.mkdtemp())
    rows = []
    for bag in sorted(bags):
        d = src.get(bag)
        if d is None:
            continue
        wins = []
        with open(d / f"{bag}.jsonl") as fh:
            for line in fh:
                wins.append(json.loads(line))
                if len(wins) >= args.sample:
                    break
        if not wins:
            continue
        alen = int(st.median(w["anchor_len"] for w in wins))
        ilen = int(st.median(ins_len[bag])) if ins_len.get(bag) else 0
        # The element minus the tnp: everything non-coding the element can hold.
        extent = max(0, ilen - alen)

        up = dn = 0
        if nc.get(bag):
            q = tmp / "nc.fa"
            q.write_text("".join(f">nc{i}\n{s}\n"
                                 for i, s in enumerate(nc[bag])))
            sj = tmp / "w.fa"
            sj.write_text("".join(
                f">w{i}\n{w['seq']}\n" for i, w in enumerate(wins)))
            r = subprocess.run(
                ["blastn", "-task", "blastn", "-word_size", "11", "-dust", "no",
                 "-evalue", "1e-5", "-max_target_seqs", "100",
                 "-query", str(q), "-subject", str(sj),
                 "-outfmt", "6 sseqid sstart send pident length"],
                capture_output=True, text=True)
            for line in r.stdout.splitlines():
                f = line.split("\t")
                if len(f) < 5 or float(f[3]) < 90:
                    continue
                w = wins[int(f[0][1:])]
                ao = w["anchor_offset"]
                s0, s1 = sorted((int(f[1]), int(f[2])))
                # anchor frame; the tnp occupies [0, anchor_len)
                a0, a1 = s0 - 1 - ao, s1 - ao
                if a1 <= 0:
                    up += 1
                elif a0 >= w["anchor_len"]:
                    dn += 1
        tot = up + dn
        fu = up / tot if tot else 0.0
        fd = dn / tot if tot else 0.0
        side_up = tot == 0 or fu >= args.min_side_frac
        side_dn = tot == 0 or fd >= args.min_side_frac
        arch = ("unknown_both" if tot == 0 else
                "[nc][tnp]+[tnp][nc]" if side_up and side_dn else
                "[nc][tnp]" if side_up else "[tnp][nc]")
        rows.append({
            "bag": bag,
            "n_sites": len(ins_len.get(bag, [])),
            "median_inserted_len": ilen,
            "median_anchor_len": alen,
            "noncoding_extent_bp": extent,
            "n_located_nc_hits": tot,
            "frac_upstream": round(fu, 3),
            "frac_downstream": round(fd, 3),
            "architecture": arch,
            "accept_upstream": side_up,
            "accept_downstream": side_dn,
            "upstream_extent_bp": (extent + args.slack) if side_up else 0,
            "downstream_extent_bp": (extent + args.slack) if side_dn else 0,
        })

    cols = list(rows[0]) if rows else []
    with open(args.out / "element_geometry.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    c = defaultdict(int)
    for r in rows:
        c[r["architecture"]] += 1
    print(f"{len(rows)} bags with geometry")
    for k, v in sorted(c.items(), key=lambda kv: -kv[1]):
        print(f"  {k:22s} {v}")
    ex = [r["noncoding_extent_bp"] for r in rows if r["noncoding_extent_bp"]]
    if ex:
        print(f"non-coding extent bp: median {st.median(ex):.0f}, "
              f"max {max(ex)}")
    print(f"-> {args.out}/element_geometry.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
