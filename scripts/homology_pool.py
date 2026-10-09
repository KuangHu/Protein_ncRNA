#!/usr/bin/env python
"""Build a remote-homology sibling pool for bags that have no family.

    homology_pool.py --targets <list> --mmseqs <m8> --phmmer <tbl> \\
                     --pfam <bag_pfam.tsv> --nr100 <bag_nr100.tsv> \\
                     --extracted <strat_windows> --out <dir>

A permutation null needs sibling clades to permute against. Eighteen bags
cleared the recalibrated floor while having none: they are singletons under the
Arm 2 family threshold (>=40% identity over >=60% of both sequences), which was
set strict on purpose so a shared domain could not merge unrelated bags. That
same strictness leaves these bags with no null at all, which is worse -- an
unchallenged score.

So the pool is rebuilt from remote homology, with two guards against making the
null easy:

  architecture   a pool member should share Pfam domain content with the target
                 where the target HAS Pfam content. For the `unknown` targets --
                 which are the most interesting ones and have no domains at all
                 -- there is nothing to match on, and that is recorded as
                 `architecture: unverifiable` rather than silently passed.
  diversity      a pool drawn from one near-identical bag is the degenerate null
                 again. Pools are required to reach a minimum number of distinct
                 source bags and a minimum distinct-protein count, and the
                 source statistics travel with every pool.

Pools are labelled `strong` / `weak` on SOURCE diversity, the same criterion
Arm 1 used, so the two are directly comparable. That label is known to be
insufficient on its own: Arm 1 called 1,497 decoy clades `strong` on source
diversity while they produced zero blocks between them. The resolving-power
label comes later, from the decoys' actual behaviour, and both are reported.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def bag(pid: str) -> str:
    return pid.split("|")[0]


def read_mmseqs(p: Path):
    with open(p) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 8:
                yield (bag(f[0]), bag(f[1]), float(f[2]), float(f[4]),
                       float(f[6]), float(f[7]))


def read_phmmer(p: Path, max_e: float) -> dict[tuple[str, str], float]:
    """phmmer --tblout: target name col 1, query name col 3, full E-value col 5.

    Returns best e-value per (query bag, hit bag). These pairs JOIN the MMseqs
    pairs rather than only annotating them: the two tools are run because they
    fail differently, and seeding candidates from one of them would make the
    other decorative.
    """
    out: dict[tuple[str, str], float] = {}
    if not p.exists():
        return out
    with open(p) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 5:
                continue
            try:
                ev = float(f[4])
            except ValueError:
                continue
            if ev > max_e:
                continue
            k = (bag(f[2]), bag(f[0]))
            if k not in out or ev < out[k]:
                out[k] = ev
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", required=True, type=Path)
    ap.add_argument("--mmseqs", required=True, type=Path)
    ap.add_argument("--phmmer", type=Path)
    ap.add_argument("--pfam", required=True, type=Path)
    ap.add_argument("--nr100", required=True, type=Path)
    ap.add_argument("--extracted", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--min-sources", type=int, default=4,
                    help="distinct source bags a pool needs; matches Arm 1's "
                         "--min-source-clades so the labels are comparable")
    ap.add_argument("--min-pool-distinct", type=int, default=200)
    ap.add_argument("--max-evalue", type=float, default=1e-3)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    targets = args.targets.read_text().split()
    nr = {r["bag_id"]: r for r in csv.DictReader(
        open(args.nr100), delimiter="\t")}
    doms = {}
    for r in csv.DictReader(open(args.pfam), delimiter="\t"):
        doms[r["bag"]] = {d.split("(")[0] for d in r["all_domains"].split(";")
                          if d}

    have_windows = set()
    if args.extracted and args.extracted.is_dir():
        for d in args.extracted.iterdir():
            if (d / f"{d.name}.fasta").exists():
                have_windows.add(d.name)

    tset = set(targets)
    ph = read_phmmer(args.phmmer, args.max_evalue) if args.phmmer else {}

    # Best hit per (target, candidate) bag pair, over all their representatives.
    best: dict[tuple[str, str], tuple[float, float, float, float]] = {}
    for q, t, pid, ev, qc, tc in read_mmseqs(args.mmseqs):
        if q == t or q not in tset or ev > args.max_evalue:
            continue
        k = (q, t)
        if k not in best or ev < best[k][1]:
            best[k] = (pid, ev, qc, tc)

    # Union, not intersection. A pair seen only by phmmer has no MMseqs
    # alignment to report, so its identity and coverage are blank rather than
    # zero -- a 0.0 there would read as a terrible alignment instead of none.
    for (q, t) in ph:
        if q in tset and q != t and (q, t) not in best:
            best[(q, t)] = (float("nan"), ph[(q, t)],
                            float("nan"), float("nan"))

    def num(x):
        return "" if x != x else x      # NaN -> blank

    rows, pools = [], {}
    for tgt in targets:
        tdom = doms.get(tgt, set())
        cand = []
        for (q, t), (pid, ev, qc, tc) in best.items():
            if q != tgt:
                continue
            shared = tdom & doms.get(t, set())
            arch = ("unverifiable" if not tdom else
                    "shared" if shared else "divergent")
            found = ("both" if (tgt, t) in ph and pid == pid else
                     "phmmer_only" if (tgt, t) in ph else "mmseqs_only")
            cand.append({
                "target": tgt, "pool_bag": t,
                "pident": num(round(pid, 1) if pid == pid else pid),
                "evalue": ev,
                "qcov": num(round(qc, 3) if qc == qc else qc),
                "tcov": num(round(tc, 3) if tc == tc else tc),
                "found_by": found,
                "phmmer_confirms": (tgt, t) in ph,
                "architecture": arch,
                "shared_domains": ";".join(sorted(shared)),
                "n_distinct": int(nr.get(t, {}).get("n_distinct", 0) or 0),
                "has_windows": t in have_windows,
            })
        # Same-architecture members first where architecture is checkable, then
        # by homology strength. A divergent-architecture member is kept but
        # ranked last and flagged -- dropping it outright would empty the pools
        # of the `unknown` targets, which have no architecture to share.
        cand.sort(key=lambda c: (
            {"shared": 0, "unverifiable": 1, "divergent": 2}[c["architecture"]],
            -int(c["phmmer_confirms"]), c["evalue"],
            -(c["pident"] if isinstance(c["pident"], float) else 0.0)))
        rows += cand

        n_src = len(cand)
        pool_distinct = sum(c["n_distinct"] for c in cand)
        n_shared = sum(1 for c in cand if c["architecture"] == "shared")
        # The Arm 1 criterion, verbatim: enough distinct sources, and no single
        # source dominating. `max_source_fraction` is computed over distinct
        # proteins, which is what a decoy is actually drawn from.
        top = max((c["n_distinct"] for c in cand), default=0)
        frac = round(top / pool_distinct, 4) if pool_distinct else 1.0
        weak = (n_src < args.min_sources
                or pool_distinct < args.min_pool_distinct
                or frac > 0.5)
        pools[tgt] = {
            "target": tgt,
            "target_distinct": int(nr.get(tgt, {}).get("n_distinct", 0) or 0),
            "target_domains": ";".join(sorted(tdom)) or "",
            "n_source_bags": n_src,
            "n_source_bags_with_windows": sum(1 for c in cand
                                              if c["has_windows"]),
            "n_shared_architecture": n_shared,
            "n_phmmer_confirmed": sum(1 for c in cand if c["phmmer_confirms"]),
            "n_found_by_one_tool_only": sum(
                1 for c in cand if c["found_by"] != "both"),
            "pool_distinct_proteins": pool_distinct,
            "max_source_fraction": frac,
            "architecture_check": ("unverifiable" if not tdom else
                                   "shared" if n_shared else "divergent_only"),
            # Source-diversity label only. Resolving power is measured after the
            # decoys are scored and is NOT implied by this.
            "null_strength_source": "weak" if weak else "strong",
        }

    cols = ["target", "pool_bag", "pident", "evalue", "qcov", "tcov",
            "found_by", "phmmer_confirms", "architecture", "shared_domains",
            "n_distinct", "has_windows"]
    with open(args.out / "pool_members.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    pcols = list(next(iter(pools.values()))) if pools else []
    with open(args.out / "pool_summary.tsv", "w") as fh:
        fh.write("\t".join(pcols) + "\n")
        for t in targets:
            fh.write("\t".join(str(pools[t][c]) for c in pcols) + "\n")

    # Pseudo-families for permutation_null: it derives the family from the text
    # before the first `__`, so each target plus its pool becomes one family.
    #
    # The relation is many-to-many and must stay that way. 46 of 136 pool bags
    # are remote homologs of more than one target -- CDS02588 serves five. A
    # bag -> family map would assign each of them to a single family and
    # silently delete it from the others, shrinking pools that were sized and
    # labelled on the full membership. So the output is (family, bag) PAIRS,
    # and the staging step writes one clade file per pair.
    pairs = []
    for i, t in enumerate(targets):
        fam = f"H{i:03d}"
        pairs.append((fam, t, "target"))
        pairs += [(fam, c["pool_bag"], "pool")
                  for c in rows if c["target"] == t]
    size = defaultdict(int)
    for f, _, _ in pairs:
        size[f] += 1
    with open(args.out / "pseudo_family_pairs.tsv", "w") as fh:
        fh.write("family\tbag_id\trole\tfamily_size\n")
        for f, b, role in sorted(pairs):
            fh.write(f"{f}\t{b}\t{role}\t{size[f]}\n")
    # Flat map kept for anything that only needs a single assignment; it is
    # lossy by construction and nothing that sizes a pool may read it.
    flat = {}
    for f, b, _ in pairs:
        flat.setdefault(b, f)
    (args.out / "pseudo_families.tsv").write_text(
        "bag_id\tfamily\tfamily_size\tn_anchors\n" + "".join(
            f"{b}\t{f}\t{size[f]}\t0\n" for b, f in sorted(flat.items())))
    needed = sorted({c["pool_bag"] for c in rows if not c["has_windows"]})
    (args.out / "pool_bags_needing_windows.txt").write_text(
        "\n".join(needed) + ("\n" if needed else ""))
    (args.out / "pool_summary.json").write_text(json.dumps(
        {"targets": pools, "n_pool_bags_needing_windows": len(needed)},
        indent=2))

    print(f"{'target':10s} {'src':>4s} {'win':>4s} {'arch':>13s} "
          f"{'distinct':>8s} {'maxfrac':>7s}  label")
    for t in targets:
        p = pools[t]
        print(f"{t:10s} {p['n_source_bags']:4d} "
              f"{p['n_source_bags_with_windows']:4d} "
              f"{p['architecture_check']:>13s} "
              f"{p['pool_distinct_proteins']:8d} "
              f"{p['max_source_fraction']:7.3f}  "
              f"{p['null_strength_source']}")
    nw = sum(1 for p in pools.values()
             if p["null_strength_source"] == "weak")
    print(f"\n{len(targets) - nw} strong / {nw} weak pools by SOURCE diversity "
          f"(Arm 1 criterion). Resolving power is not implied by this label and "
          f"is measured after the decoys are scored.")
    print(f"{len(needed)} pool bags still need window extraction")
    print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
