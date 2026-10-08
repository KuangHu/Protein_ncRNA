#!/usr/bin/env python
"""Step 3: clade decomposition and per-clade feasibility.

A family-wide alignment of 138,897 IS110 transposases is biological soup, and a
conserved block found in it says nothing. Step 4 needs *alignable subfamilies*:
tight enough that a shared RNA is really shared, loose enough that substitutions
exist for covariation to score.

    protein identity   expected behaviour
    > 95%              RNA nearly identical; R-scape has no substitutions
    70-90%             usually good for a tight element subfamily
    35-70%             may keep the scaffold, may lose alignability
    < 35%              normally too remote unless the RNA is severely constrained

The sweet spot differs per family, so this runs a ladder (90/70/50/35) and lets
the selection rule pick per family rather than imposing one threshold.

    clade_decompose.py --census <corrected_census> --out out/clades

Outputs `cluster_members.tsv`, `clade_summary.tsv` and `selected_clades.tsv`.
The last is the entry point for step 4.
"""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_ncrna.compute import require_slurm
from protein_ncrna.diversity import hash_sample
from protein_ncrna.seqio import read_fasta
from protein_ncrna.tools import which

import contextlib


@contextlib.contextmanager
def _null_writer():
    """Stand-in for the cluster_members handle when clusters are being reused."""
    class _N:
        def write(self, *_a):
            pass
    yield _N()


LADDER = [0.90, 0.70, 0.50, 0.35]

# Families to decompose, in the priority order the census results justify.
# TnpB first: omegaRNA is a known, published guide RNA encoded immediately
# downstream of tnpB, so anchor-to-RNA geometry is the clearest of any family.
DEFAULT_FAMILIES = ["TnpB_Cas12f", "IS110", "GroupII_RT", "RT", "Cas3"]

# A clade is worth taking to step 4 when it is neither a strain pileup nor a
# grab-bag. Bounds are deliberately wide; the per-clade table is the real output.
SELECT = {
    # Required.
    "min_nr100": 8,
    "min_complete_windows": 8,
    "min_median_identity": 0.35,
    "max_median_identity": 0.90,
    # Soft sanity only: stops a clade that is almost entirely contig edge. What
    # step 4 needs is how many complete +/-5 kb windows can be aligned, which is
    # the count above, not the proportion. Requiring frac >= 0.8 discarded ~70%
    # of good IS110 clades, because mobile elements cause the assembly breaks
    # that truncate their own windows.
    "min_complete_window_frac": 0.50,
    # Reported, not enforced: clade size and species breadth are recorded so an
    # oversized or single-species clade is visible, but neither rejects on its
    # own. The identity window already bounds alignability.
    "report_max_nr100": 200,
    "report_min_species": 2,
}


def run_mmseqs_cluster(faa: Path, ident: float, tmp: Path, threads: int) -> dict[str, str]:
    """Cluster one family FASTA. Returns anchor_id -> cluster representative id."""
    mmseqs = which("mmseqs")
    pref = tmp / f"c{int(ident * 100)}"
    pref.parent.mkdir(parents=True, exist_ok=True)
    cmd = [mmseqs, "easy-cluster", str(faa), str(pref), str(tmp / f"t{int(ident * 100)}"),
           "--min-seq-id", str(ident), "-c", "0.8", "--cov-mode", "0",
           "--threads", str(threads), "-v", "1"]
    proc = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        raise RuntimeError(f"mmseqs failed on {faa.name} at {ident}: "
                           f"{proc.stderr.decode(errors='replace')[-300:]}")
    members: dict[str, str] = {}
    with open(str(pref) + "_cluster.tsv") as fh:
        for line in fh:
            rep, mem = line.rstrip("\n").split("\t")[:2]
            members[mem] = rep
    return members


def approx_median_identity(faa: Path, members: list[str], seqs: dict[str, str],
                           tmp: Path, threads: int, cap: int = 300,
                           seed: int = 1) -> float | None:
    """Median pairwise AA identity within a clade, from an all-vs-all search.

    Capped: identity is a property of the clade, and a few hundred sequences
    estimate it as well as thousands at a fraction of the cost.

    The sample is a seeded hash of the member ids, not the first `cap` in MMseqs
    output order. MMseqs emits members grouped by representative and ordered by
    the input FASTA, which tracks accession and submission batch, so the head of
    the list is not a neutral slice of the clade.

    It is deliberately *not* a diversity-spread sample, unlike the window cap.
    This number gates clade selection as an estimate of how similar the clade's
    members typically are; a farthest-point subset is enriched for outliers by
    construction and would bias the median down, increasingly so for the large
    clades where the cap actually binds. Spread the sample when the subset must
    represent the clade's breadth, and hash it when the subset must stand in for
    the clade's typical member.

    The sample is drawn from the clade's *distinct* sequences. Without that it
    is not an estimate of anything: these clades are mostly duplicates -- one
    GroupII_RT clade has 799 members and 29 distinct proteins, a ratio of 0.036
    -- so a sample of raw members is overwhelmingly duplicate-vs-duplicate pairs
    and the median converges on 1.0 regardless of how diverged the clade really
    is. That is what it did: the same clade measured 0.451 from an ordered
    sample and 1.000 from a random one, and the 1.000 pushed it past the
    `identity > 0.9` gate, dropping GroupII_RT from the 0.35 rung to the 0.70
    rung and 3,415 alignable windows to 217. Neither number was a measurement.
    """
    if len(members) < 2:
        return None
    # Deduplicate on sequence, keeping the first id for each, so the estimate
    # measures divergence rather than copy number.
    first: dict[str, str] = {}
    for m in sorted(members):
        if m in seqs:
            first.setdefault(seqs[m], m)
    uniq = sorted(first.values())
    if len(uniq) < 2:
        # A clade of one distinct sequence is perfectly conserved, and saying so
        # is more useful than returning None and skipping the gate.
        return 1.0
    sample = hash_sample(uniq, cap, seed)
    tmp.mkdir(parents=True, exist_ok=True)
    q = tmp / "clade.faa"
    with open(q, "w") as fh:
        for m in sample:
            fh.write(f">{m}\n{seqs[m]}\n")
    res = tmp / "hits.m8"
    proc = subprocess.run(
        [which("mmseqs"), "easy-search", str(q), str(q), str(res), str(tmp / "t"),
         "--threads", str(threads), "-v", "1", "--max-seqs", "400",
         "--format-output", "query,target,fident"],
        stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    if proc.returncode != 0 or not res.exists():
        return None
    vals = []
    with open(res) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[0] != f[1]:
                vals.append(float(f[2]))
    return round(statistics.median(vals), 4) if vals else None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--families", nargs="*", default=DEFAULT_FAMILIES)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--identity", action="store_true",
                    help="also estimate median pairwise identity per selected clade")
    ap.add_argument("--reuse-clusters", action="store_true",
                    help="reuse an existing cluster_members.tsv instead of "
                         "re-running mmseqs; for iterating on the selection rule")
    ap.add_argument("--neighborhood", type=Path, default=None,
                    help="neighborhood.tsv from neighborhood_scan.py, for Cas3 "
                         "verification and RT subtyping")
    args = ap.parse_args()
    require_slurm("clade decomposition")
    args.out.mkdir(parents=True, exist_ok=True)

    # Per-anchor metadata from the corrected census.
    meta: dict[str, dict] = {}
    with open(args.census / "anchors.tsv") as fh:
        cols = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            meta[f[cols["anchor_id"]]] = {
                "family": f[cols["family"]],
                "species": f[cols["species"]],
                "genome": f[cols["genome"]],
                "model": f[cols["model"]],
                "window_complete": f[cols["window_complete"]] == "True",
                "anchor_partial": f[cols["anchor_partial"]] == "True",
                "truncated": (f[cols["truncated_left"]] == "True"
                              or f[cols["truncated_right"]] == "True"),
                "aa_len": int(f[cols["anchor_len"]]) // 3,
            }
    print(f"loaded metadata for {len(meta)} anchors")

    nbr = load_neighborhood(args.neighborhood) if args.neighborhood else {}
    if nbr:
        print(f"loaded neighborhood calls for {len(nbr)} anchors")

    tmp_root = args.out / "_tmp"
    cm_path = args.out / "cluster_members.tsv"
    rows_summary = []
    seqs_by_fam: dict[str, dict[str, str]] = {}

    # Pairwise identity is the expensive column; carry it across re-runs so the
    # selection rule can be iterated without paying for it again.
    prior_identity: dict[tuple, float] = {}
    cs_path = args.out / "clade_summary.tsv"
    if args.reuse_clusters and cs_path.exists():
        with open(cs_path) as fh:
            c = {k: i for i, k in enumerate(next(fh).rstrip("\n").split("\t"))}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                v = f[c["median_identity"]] if "median_identity" in c else ""
                if v:
                    prior_identity[(f[c["family"]], float(f[c["threshold"]]),
                                    f[c["cluster_id"]])] = float(v)
        print(f"recovered {len(prior_identity)} identity estimates")

    if args.reuse_clusters and cm_path.exists():
        print(f"reusing clusters from {cm_path}")
        for fam in args.families:
            faa = args.census / "proteins" / f"{fam}.faa"
            if faa.exists():
                seqs_by_fam[fam] = {h.split()[0]: s for h, s in read_fasta(str(faa))}
        grouped: dict[tuple, list[str]] = defaultdict(list)
        with open(cm_path) as cm:
            next(cm)
            for line in cm:
                f = line.rstrip("\n").split("\t")
                grouped[(f[0], float(f[1]), f[2])].append(f[3])
        for (fam, ident, rep), mem in grouped.items():
            if fam in seqs_by_fam:
                r = summarize_clade(fam, ident, rep, mem, meta, seqs_by_fam[fam], nbr)
                r["median_identity"] = prior_identity.get((fam, ident, rep))
                rows_summary.append(r)
        print(f"  {len(rows_summary)} clades reloaded")

    with open(cm_path, "w") if not args.reuse_clusters else _null_writer() as cm:
        if not args.reuse_clusters:
            cm.write("family\tthreshold\tcluster_id\tanchor_id\n")
        for fam in ([] if args.reuse_clusters else args.families):
            faa = args.census / "proteins" / f"{fam}.faa"
            if not faa.exists():
                print(f"  {fam}: no FASTA, skipping", file=sys.stderr)
                continue
            seqs = {h.split()[0]: s for h, s in read_fasta(str(faa))}
            seqs_by_fam[fam] = seqs
            print(f"\n[{fam}] {len(seqs)} anchors")

            for ident in LADDER:
                members = run_mmseqs_cluster(faa, ident, tmp_root / fam, args.threads)
                clusters: dict[str, list[str]] = defaultdict(list)
                for aid, rep in members.items():
                    clusters[rep].append(aid)
                    cm.write(f"{fam}\t{ident}\t{rep}\t{aid}\n")
                print(f"  {int(ident*100)}%: {len(clusters)} clusters")

                for rep, mem in clusters.items():
                    rows_summary.append(summarize_clade(fam, ident, rep, mem, meta,
                                                        seqs, nbr))

    # Per-clade feasibility table.
    cols_out = ["family", "threshold", "cluster_id", "n_anchors", "n_nr100", "n_species",
                "n_genomes", "median_aa_len", "n_complete_windows",
                "complete_window_frac", "contig_edge_frac",
                "duplicate_frac", "rt_subtype", "cas_neighborhood", "median_identity",
                "flags", "selected", "reject_reason"]

    # Apply the cheap criteria first so that the expensive identity estimate is
    # paid only for clades that could still be selected. Without this pass
    # reject_reason is still empty and every clade looks like a candidate.
    for r in rows_summary:
        apply_selection(r)

    if args.identity:
        cand = [r for r in rows_summary if r["selected"]]
        print(f"\nestimating pairwise identity for {len(cand)} candidate clades")
        by_key = {(r["family"], r["threshold"], r["cluster_id"]): r for r in rows_summary}
        with open(cm_path) as fh:
            next(fh)
            mem_of: dict[tuple, list[str]] = defaultdict(list)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                key = (f[0], float(f[1]), f[2])
                if key in by_key and by_key[key]["selected"]:
                    mem_of[key].append(f[3])
        for i, (key, mem) in enumerate(mem_of.items(), 1):
            r = by_key[key]
            r["median_identity"] = approx_median_identity(
                None, mem, seqs_by_fam[key[0]], tmp_root / "ident", args.threads)
            if i % 50 == 0:
                print(f"  {i}/{len(mem_of)}", flush=True)

    for r in rows_summary:
        apply_selection(r)

    with open(args.out / "clade_summary.tsv", "w") as fh:
        fh.write("\t".join(cols_out) + "\n")
        for r in sorted(rows_summary, key=lambda r: (r["family"], -r["threshold"],
                                                     -r["n_anchors"])):
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c, ""))
                               for c in cols_out) + "\n")

    sel = [r for r in rows_summary if r["selected"]]

    def _write(path: Path, rs: list[dict]) -> None:
        with open(path, "w") as fh:
            fh.write("\t".join(cols_out) + "\n")
            for r in sorted(rs, key=lambda r: (r["family"], -r["n_nr100"])):
                fh.write("\t".join("" if r.get(c) is None else str(r.get(c, ""))
                                   for c in cols_out) + "\n")

    # Clades at different rungs are nested: the same locus is selected two or
    # three times over. Step 4 must consume exactly one rung per family or it
    # aligns the same block repeatedly and reports the duplicates as support.
    rec = recommend_thresholds(rows_summary)
    (args.out / "recommended_thresholds.json").write_text(json.dumps(rec, indent=2))
    final = [r for r in sel if r["threshold"] == rec[r["family"]]["recommended_threshold"]]

    _write(args.out / "selected_clades_all_thresholds.tsv", sel)
    _write(args.out / "selected_clades.tsv", final)

    print(f"\n{'family':<14}{'rung':>6}{'clades':>8}{'nr100':>8}{'windows':>10}{'identity':>10}")
    for fam, v in sorted(rec.items()):
        print(f"{fam:<14}{v['recommended_threshold']:>6}{v['n_selected']:>8}"
              f"{v['nr100_in_selected']:>8}{v['complete_windows_in_selected']:>10}"
              f"{str(v['median_identity_of_selected']):>10}")
    print(f"\n{len(final)} clades at the recommended rung "
          f"(of {len(sel)} selected across all rungs)")

    report_summary(rows_summary, final, args.out)
    return 0


def summarize_clade(fam: str, ident: float, rep: str, mem: list[str],
                    meta: dict, seqs: dict, nbr: dict) -> dict:
    m = [meta[a] for a in mem if a in meta]
    if not m:
        m = [{"species": "unknown", "genome": "?", "window_complete": False,
              "truncated": True, "aa_len": 0, "model": "?"}]
    uniq_seq = len({seqs[a] for a in mem if a in seqs})
    genomes = Counter(x["genome"] for x in m)
    r = {
        "family": fam,
        "threshold": ident,
        "cluster_id": rep,
        "n_anchors": len(mem),
        "n_nr100": uniq_seq,
        "n_species": len({x["species"] for x in m}),
        "n_genomes": len(genomes),
        "median_aa_len": int(statistics.median([x["aa_len"] for x in m])),
        "n_complete_windows": sum(x["window_complete"] for x in m),
        "complete_window_frac": round(sum(x["window_complete"] for x in m) / len(m), 4),
        "contig_edge_frac": round(sum(x["truncated"] for x in m) / len(m), 4),
        # Many copies per genome is normal for a mobile element, but an extreme
        # value flags a repeat family rather than a discrete locus.
        "duplicate_frac": round(1 - len(genomes) / len(mem), 4),
        "median_identity": None,
        "rt_subtype": "",
        "cas_neighborhood": "",
        "reject_reason": "",
        "flags": "",
        "selected": False,
    }
    if fam in ("RT", "GroupII_RT"):
        r["rt_subtype"] = rt_subtype(fam, m, mem, nbr)
    if fam == "Cas3":
        r["cas_neighborhood"] = cas_neighborhood_call(mem, nbr)
    return r


def rt_subtype(fam: str, m: list[dict], mem: list[str], nbr: dict) -> str:
    """Coarse RT labels, only enough to stop incompatible loci being pooled.

    Deliberately not a retron caller: 'retron' is not claimed from RT homology
    alone, since generic RT covers retrons, Abi-like defence RTs, DGRs, group II
    remnants and host enzymes alike.
    """
    if fam == "GroupII_RT":
        return "groupII"
    med_len = statistics.median([x["aa_len"] for x in m])
    hits = Counter()
    for a in mem:
        for g in nbr.get(a, ()):
            hits[g] += 1
    frac = {g: n / len(mem) for g, n in hits.items()}
    defence = {"HNH", "TOPRIM", "HEPN", "Abi", "SIR2", "TIR"}
    if any(frac.get(g, 0) > 0.5 for g in defence):
        return "candidate_defence"
    if med_len >= 600:
        return "long_maturase_like"
    if 250 <= med_len <= 450:
        return "unknown_compact"
    return "unknown"


def cas_neighborhood_call(mem: list[str], nbr: dict) -> str:
    """Cas3 is provisional until cascade genes are seen beside it.

    PF18019 is an HD domain and the HD clan is promiscuous, so an HD-only
    singleton is downweighted rather than trusted.
    """
    if not nbr:
        return "unscanned"
    cas = {"Cas5", "Cas6", "Cas7", "Cas1", "Cas2", "Cascade", "Cas8", "Cas11"}
    with_cas = sum(1 for a in mem if cas & set(nbr.get(a, ())))
    frac = with_cas / len(mem)
    if frac >= 0.5:
        return f"verified ({frac:.2f})"
    if frac >= 0.1:
        return f"partial ({frac:.2f})"
    return f"hd_only ({frac:.2f})"


def apply_selection(r: dict) -> None:
    """Flag a clade as a step-4 entry point, recording why when it is not."""
    why = []
    if r["n_nr100"] < SELECT["min_nr100"]:
        why.append(f"nr100<{SELECT['min_nr100']}")
    if r["n_complete_windows"] < SELECT["min_complete_windows"]:
        why.append(f"complete_windows<{SELECT['min_complete_windows']}")
    if r["complete_window_frac"] < SELECT["min_complete_window_frac"]:
        why.append(f"complete_window_frac<{SELECT['min_complete_window_frac']}")
    mi = r.get("median_identity")
    if mi is not None:
        if mi < SELECT["min_median_identity"]:
            why.append(f"identity<{SELECT['min_median_identity']}")
        if mi > SELECT["max_median_identity"]:
            why.append(f"identity>{SELECT['max_median_identity']}")
    # Cas3 is a *provisional* anchor family: PF18019 is an HD domain and the HD
    # clan is promiscuous, so a Cas3 clade is not an entry point to step 4 until
    # cascade genes are actually seen beside it. Unscanned is not a pass.
    if r["family"] == "Cas3":
        call = r["cas_neighborhood"]
        if call.startswith("hd_only"):
            why.append("cas3 HD-only, no cascade neighbours")
        elif call in ("", "unscanned"):
            why.append("cas3 unverified: run neighborhood_scan.py first")
    r["reject_reason"] = "; ".join(why)
    r["selected"] = not why
    # Advisory only — these never reject, but step 4 should see them.
    flags = []
    if r["n_nr100"] > SELECT["report_max_nr100"]:
        flags.append(f"large_clade(nr100={r['n_nr100']})")
    if r["n_species"] < SELECT["report_min_species"]:
        flags.append("single_species")
    r["flags"] = ",".join(flags)


def load_neighborhood(path: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    if not path.exists():
        print(f"neighborhood file {path} not found", file=sys.stderr)
        return {}
    with open(path) as fh:
        cols = {c: i for i, c in enumerate(next(fh).rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[cols["anchor_id"]]].add(f[cols["neighbor_group"]])
    return out


def recommend_thresholds(rows: list[dict]) -> dict[str, dict]:
    """Pick one ladder rung per family.

    The ladder exists because the right identity cut differs per family — TnpB
    may be alignable at 50-70% while IS110 needs tighter because its RNAs churn.
    Clades at different rungs are nested, so step 4 must consume exactly one
    rung per family or it will process the same locus repeatedly.

    Chosen by the most in-spec clades. Ties break toward the rung whose selected
    clades are *least* identical, because both rungs already passed the same
    0.35-0.90 band and within that band the only thing that separates them is
    covariation power. Breaking ties toward the tighter cut instead would, for
    GroupII_RT here, trade a 0.59-identity rung for a 0.79-identity one and drop
    3415 alignable windows to 217 — the silent-failure direction.
    """
    by_fam: dict[str, dict[float, list]] = defaultdict(lambda: defaultdict(list))
    for r in rows:
        by_fam[r["family"]][r["threshold"]].append(r)

    out = {}
    for fam, rungs in by_fam.items():
        scored = []
        for thr, rs in rungs.items():
            sel = [r for r in rs if r["selected"]]
            idents = [r["median_identity"] for r in sel
                      if r.get("median_identity") is not None]
            scored.append({
                "threshold": thr,
                "n_selected": len(sel),
                "n_clades": len(rs),
                "nr100_in_selected": sum(r["n_nr100"] for r in sel),
                "complete_windows_in_selected": sum(r["n_complete_windows"]
                                                    for r in sel),
                "median_identity_of_selected": (round(statistics.median(idents), 4)
                                                if idents else None),
                "median_species": (statistics.median([r["n_species"] for r in sel])
                                   if sel else 0),
            })
        best = max(scored, key=lambda s: (s["n_selected"],
                                          -(s["median_identity_of_selected"] or 1.0)))
        out[fam] = {"recommended_threshold": best["threshold"],
                    "n_selected": best["n_selected"],
                    "nr100_in_selected": best["nr100_in_selected"],
                    "complete_windows_in_selected": best["complete_windows_in_selected"],
                    "median_identity_of_selected": best["median_identity_of_selected"],
                    "ladder": sorted(scored, key=lambda s: -s["threshold"])}
    return out


def report_summary(rows: list[dict], sel: list[dict], out: Path) -> None:
    print(f"\n{'family':<14}{'thr':>6}{'clades':>8}{'selected':>10}{'anchors in selected':>21}")
    agg: dict[tuple, dict] = defaultdict(lambda: {"n": 0, "sel": 0, "anch": 0})
    for r in rows:
        a = agg[(r["family"], r["threshold"])]
        a["n"] += 1
        a["sel"] += r["selected"]
        a["anch"] += r["n_anchors"] if r["selected"] else 0
    for (fam, thr), a in sorted(agg.items()):
        print(f"{fam:<14}{thr:>6}{a['n']:>8}{a['sel']:>10}{a['anch']:>21}")

    reasons = Counter(r["reject_reason"].split(";")[0].strip()
                      for r in rows if r["reject_reason"])
    print("\ntop rejection reasons:")
    for k, v in reasons.most_common(8):
        print(f"  {v:>7}  {k}")

    (out / "selection_rule.json").write_text(json.dumps(
        {"ladder": LADDER, "rule": SELECT, "n_clades": len(rows),
         "n_selected": len(sel)}, indent=2))
    print(f"\n{len(sel)} clades at the recommended rung -> {out}/selected_clades.tsv")


if __name__ == "__main__":
    raise SystemExit(main())
