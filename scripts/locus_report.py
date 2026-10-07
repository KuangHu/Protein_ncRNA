#!/usr/bin/env python
"""Per-locus report for one or more step-4 blocks.

    locus_report.py --blocks <k3_dir> --windows <windows_dir> --out <dir> \\
                    --block-id <BLOCK> [--block-id <BLOCK> ...] \\
                    [--fold <fold_dir>]

Step 4 says a block recurs at a tight relative coordinate; it does not say what
the block *is*. This report lays out the locus around each member so that can be
argued from the data: where the block sits relative to the anchor ORF, which
predicted genes flank it, how far the surrounding noncoding run extends, and
whether the clade's other blocks are neighbours of this one or unrelated.

That last column is the one that separates the interpretations for a group II
RT. Two blocks a few hundred bp apart on the same side of the anchor, inside one
contiguous noncoding run, are two arms of a single element; two blocks on
opposite sides are not. `run_start`/`run_end` give the run each block lives in,
and `other_blocks_in_run` names the company it keeps.

Blind: ORF geometry and the block's own members only. No curated RNA family,
no intron model, no Rfam. Whether a block *is* a group II intron arm is a
question for the benchmark branch, not for this file -- what is inferable here
is geometry, and geometry is reported as geometry.

Outputs per block: `<id>.locus_table.tsv`, `<id>.relative_coordinate_plot.svg`,
`<id>.summary.json`, and, when `--fold` points at a fold_blocks.py run,
`<id>.alignment.sto` and `<id>.consensus_structure.txt` copied beside them.
"""

from __future__ import annotations

import argparse
import json
import shutil
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


def write_tsv(path: Path, rows: list[dict], cols: list[str]) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def noncoding_runs(mask: list[tuple[int, int]], win_len: int) -> list[tuple[int, int]]:
    """Complement of the merged coding mask, in window coordinates."""
    out, prev = [], 0
    for ms, me in sorted(mask):
        if ms > prev:
            out.append((prev, ms))
        prev = max(prev, me)
    if prev < win_len:
        out.append((prev, win_len))
    return out


def flanks(s: int, e: int, mask: list[tuple[int, int]]) -> dict:
    """Nearest predicted ORF on each side, in window coordinates."""
    up = [(ms, me) for ms, me in mask if me <= s]
    dn = [(ms, me) for ms, me in mask if ms >= e]
    u = max(up, key=lambda t: t[1]) if up else None
    d = min(dn, key=lambda t: t[0]) if dn else None
    return {"upstream_orf_end": u[1] if u else "", "dist_upstream_orf": s - u[1] if u else "",
            "downstream_orf_start": d[0] if d else "", "dist_downstream_orf": d[0] - e if d else ""}


def pid(a: str, b: str) -> float:
    n = m = 0
    for x, y in zip(a, b):
        if x == "-" and y == "-":
            continue
        n += 1
        m += x.upper() == y.upper()
    return m / n if n else 0.0


def svg_plot(path: Path, rows: list[dict], bid: str, others: list[dict]) -> None:
    """Members as horizontal bars in anchor-frame coordinates, anchor ORF at x=0.

    Hand-rolled SVG: the plot is a few hundred rectangles on one axis, and a
    matplotlib dependency on the cluster side is not worth that.
    """
    rows = sorted(rows, key=lambda r: r["relative_start"])
    lo = min([r["run_start_rel"] for r in rows] + [r["relative_start"] for r in rows] + [0])
    hi = max([r["run_end_rel"] for r in rows] + [r["relative_end"] for r in rows]
             + [rows[0]["anchor_len"]])
    pad = max(50, (hi - lo) // 20)
    lo, hi = lo - pad, hi + pad
    W, left, top = 1000.0, 70.0, 70.0
    rh = max(1.0, min(6.0, 520.0 / max(1, len(rows))))
    H = top + rh * len(rows) + 60
    sx = lambda x: left + (x - lo) / (hi - lo) * (W - left - 20)

    p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" height="{H:.0f}" '
         f'font-family="sans-serif" font-size="11">',
         f'<rect width="{W:.0f}" height="{H:.0f}" fill="white"/>',
         f'<text x="{left}" y="22" font-size="14">{bid}</text>',
         f'<text x="{left}" y="40" fill="#666">{len(rows)} members &#183; '
         f'anchor-frame coordinates, x=0 is the anchor start codon</text>']
    # anchor ORF band
    p.append(f'<rect x="{sx(0):.1f}" y="{top - 14:.1f}" '
             f'width="{sx(rows[0]["anchor_len"]) - sx(0):.1f}" height="{rh*len(rows)+14:.1f}" '
             f'fill="#e8e8e8"/>')
    p.append(f'<text x="{sx(0):.1f}" y="{top - 20:.1f}" fill="#777">anchor ORF</text>')
    for o in others:
        if o["block_id"] == bid:
            continue
        p.append(f'<rect x="{sx(o["start"]):.1f}" y="{top - 8:.1f}" '
                 f'width="{max(1.0, sx(o["end"]) - sx(o["start"])):.1f}" '
                 f'height="{rh*len(rows)+8:.1f}" fill="#f2c8c8" opacity="0.5"/>')
        p.append(f'<text x="{sx(o["start"]):.1f}" y="{H-32:.1f}" fill="#b44" '
                 f'font-size="9">{o["block_id"].rsplit("__", 1)[-1]}</text>')
    for i, r in enumerate(rows):
        y = top + i * rh
        p.append(f'<rect x="{sx(r["run_start_rel"]):.1f}" y="{y + rh*0.35:.2f}" '
                 f'width="{max(0.5, sx(r["run_end_rel"]) - sx(r["run_start_rel"])):.1f}" '
                 f'height="{max(0.4, rh*0.3):.2f}" fill="#cfd8e8"/>')
        c = "#2a6" if r["strand"] == "+" else "#b60"
        p.append(f'<rect x="{sx(r["relative_start"]):.1f}" y="{y:.2f}" '
                 f'width="{max(0.8, sx(r["relative_end"]) - sx(r["relative_start"])):.1f}" '
                 f'height="{rh*0.85:.2f}" fill="{c}"/>')
    for t in range(5):
        x = lo + (hi - lo) * t / 4
        p.append(f'<line x1="{sx(x):.1f}" y1="{top - 14:.1f}" x2="{sx(x):.1f}" '
                 f'y2="{top + rh*len(rows):.1f}" stroke="#ddd"/>')
        p.append(f'<text x="{sx(x):.1f}" y="{H-14:.1f}" text-anchor="middle" '
                 f'fill="#444">{int(x)}</text>')
    p.append(f'<text x="{left}" y="{H-2:.1f}" fill="#666" font-size="9">'
             f'green/orange = block on + / - strand, pale blue = contiguous '
             f'noncoding run, pink = other blocks of this clade</text></svg>')
    path.write_text("\n".join(p))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--blocks", required=True, type=Path, help="step-4 output dir")
    ap.add_argument("--windows", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--block-id", action="append", required=True)
    ap.add_argument("--fold", type=Path, help="fold_blocks.py out dir, for sto/structure")
    ap.add_argument("--tmp", type=Path, default=None)
    args = ap.parse_args()
    require_slurm("locus report")
    args.out.mkdir(parents=True, exist_ok=True)

    blocks = {b["block_id"]: b for b in read_tsv(args.blocks / "blocks.tsv")}
    members = defaultdict(list)
    for m in read_tsv(args.blocks / "block_members.tsv"):
        members[m["block_id"]].append(m)

    wanted = [b for b in args.block_id]
    missing = [b for b in wanted if b not in blocks]
    if missing:
        raise SystemExit(f"unknown block_id: {missing}")

    tmp_root = Path(tempfile.mkdtemp(dir=args.tmp))
    by_clade = defaultdict(list)
    for b in wanted:
        by_clade[blocks[b]["clade"]].append(b)

    summary = {}
    for clade, bids in by_clade.items():
        recs = _db.load_clade(args.windows / f"{clade}.fasta")
        by_id = {r["anchor_id"]: r for r in recs}
        mask = _db.call_orfs_mask(recs, tmp_root / clade[:60])
        runs = {i: noncoding_runs(mask.get(i, []), by_id[i]["win_len"]) for i in by_id}
        # Every block this clade produced, in anchor frame, for the neighbour column.
        clade_blocks = [{"block_id": k, "start": int(v["median_start"]),
                         "end": int(v["median_start"]) + int(v["median_len"])}
                        for k, v in blocks.items() if v["clade"] == clade]

        for bid in bids:
            rows = []
            for m in members[bid]:
                r = by_id.get(m["anchor_id"])
                if not r:
                    continue
                off, al = r["anchor_offset"], r["anchor_len"]
                s, e = int(m["start"]) + off, int(m["end"]) + off
                run = next(((a, b) for a, b in runs[m["anchor_id"]]
                            if a <= s and b >= e), (s, e))
                fl = flanks(s, e, mask.get(m["anchor_id"], []))
                # Signed distance to the anchor ORF: <0 upstream of the start
                # codon, >0 past the stop, 0 if the block abuts or overlaps.
                d = (int(m["end"]) if int(m["end"]) <= 0 else
                     (int(m["start"]) - al if int(m["start"]) >= al else 0))
                rows.append({
                    "block_id": bid, "anchor_id": m["anchor_id"],
                    "relative_start": int(m["start"]), "relative_end": int(m["end"]),
                    "block_len": int(m["end"]) - int(m["start"]),
                    "strand": m.get("strand", ""), "pident_to_reference": m.get("pident", ""),
                    "dist_to_anchor_orf": d,
                    "anchor_len": al, "win_len": r["win_len"],
                    "run_start_rel": run[0] - off, "run_end_rel": run[1] - off,
                    "run_len": run[1] - run[0],
                    "block_frac_of_run": round((e - s) / max(1, run[1] - run[0]), 3),
                    "upstream_orf_end_rel": (fl["upstream_orf_end"] - off
                                             if fl["upstream_orf_end"] != "" else ""),
                    "dist_upstream_orf": fl["dist_upstream_orf"],
                    "downstream_orf_start_rel": (fl["downstream_orf_start"] - off
                                                 if fl["downstream_orf_start"] != "" else ""),
                    "dist_downstream_orf": fl["dist_downstream_orf"],
                    "other_blocks_in_run": ",".join(
                        o["block_id"].rsplit("__", 1)[-1] for o in clade_blocks
                        if o["block_id"] != bid
                        and o["start"] >= run[0] - off and o["end"] <= run[1] - off)
                        or "-",
                })
            if not rows:
                print(f"  {bid}: no members resolved", file=sys.stderr)
                continue

            cols = ["block_id", "anchor_id", "relative_start", "relative_end",
                    "block_len", "strand", "pident_to_reference", "dist_to_anchor_orf",
                    "run_start_rel", "run_end_rel", "run_len", "block_frac_of_run",
                    "upstream_orf_end_rel", "dist_upstream_orf",
                    "downstream_orf_start_rel", "dist_downstream_orf",
                    "other_blocks_in_run", "anchor_len", "win_len"]
            short = bid.rsplit("__", 1)[-1]
            stem = args.out / f"{short}"
            write_tsv(Path(f"{stem}.locus_table.tsv"), rows, cols)
            svg_plot(Path(f"{stem}.relative_coordinate_plot.svg"), rows, bid, clade_blocks)

            # Pairwise identity distribution, read off the fold alignment if one
            # exists so the number matches what R-scape was given.
            ids = []
            if args.fold:
                aln = args.fold / bid / "aln.fa"
                if aln.exists():
                    from protein_ncrna.seqio import read_fasta
                    seqs = [s for _, s in read_fasta(str(aln))]
                    step = max(1, len(seqs) // 60)
                    sub = seqs[::step]
                    ids = [pid(sub[i], sub[j]) for i in range(len(sub))
                           for j in range(i + 1, len(sub))]
                for src, dst in (("aln.sto", f"{short}.alignment.sto"),
                                 ("structure.txt", f"{short}.consensus_structure.txt"),
                                 ("aln.names.tsv", f"{short}.alignment.names.tsv")):
                    p = args.fold / bid / src
                    if p.exists():
                        shutil.copy(p, args.out / dst)

            L = [r["block_len"] for r in rows]
            S = [r["relative_start"] for r in rows]
            q = lambda v, f: round(sorted(v)[min(len(v) - 1, int(f * len(v)))], 4)
            summary[bid] = {
                "family": blocks[bid]["family"], "clade": clade,
                "score": float(blocks[bid]["composite_score"]),
                "n_members": len(rows),
                "relative_start_median": int(statistics.median(S)),
                "relative_start_iqr": [q(S, .25), q(S, .75)],
                "block_len_median": int(statistics.median(L)),
                "block_len_range": [min(L), max(L)],
                "strand": dict(sorted(
                    {s: sum(1 for r in rows if r["strand"] == s)
                     for s in {r["strand"] for r in rows}}.items())),
                "side_of_anchor": {
                    "upstream": sum(1 for r in rows if r["relative_end"] <= 0),
                    "overlapping": sum(1 for r in rows if r["relative_start"] < rows[0]["anchor_len"]
                                       and r["relative_end"] > 0),
                    "downstream": sum(1 for r in rows if r["relative_start"] >= r["anchor_len"]),
                },
                "block_frac_of_noncoding_run_median":
                    round(statistics.median([r["block_frac_of_run"] for r in rows]), 3),
                "noncoding_run_len_median":
                    int(statistics.median([r["run_len"] for r in rows])),
                "other_blocks_in_same_run": sorted(
                    {x for r in rows for x in r["other_blocks_in_run"].split(",")
                     if x != "-"}),
                "pairwise_identity": ({
                    "n_pairs": len(ids), "median": round(statistics.median(ids), 4),
                    "p05": q(ids, .05), "p95": q(ids, .95),
                    "min": round(min(ids), 4), "max": round(max(ids), 4)} if ids else None),
            }
            (args.out / f"{short}.summary.json").write_text(
                json.dumps(summary[bid], indent=2))
            print(f"{bid}\n  " + "\n  ".join(
                f"{k}: {v}" for k, v in summary[bid].items()), flush=True)

    (args.out / "locus_report.json").write_text(json.dumps(summary, indent=2))
    print(f"\n-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
