#!/usr/bin/env python
"""Step 6 — comparative RNA structure test on step-4 candidate blocks.

    fold_blocks.py --fasta <export_dir> --out <dir> [--min-score 0.90]
                   [--exclude-benchmark-only]

MAFFT -> RNAalifold -> R-scape, on the padded FASTA of each block. Nothing here
reads RetronDB, Rfam, msr/msd or any curated covariance model: a block is judged
by whether *its own* members covary, not by resemblance to a known family.

What the numbers mean, and do not
---------------------------------
RNAalifold will return a structure for anything, including a random alignment.
The consensus MFE on its own is close to meaningless; what matters is whether
compensatory substitutions support the helices, which is what R-scape tests.

R-scape needs diverged sequences. An alignment of near-identical members yields
"no significant covariation" because there were no substitutions to covary, not
because there is no structure -- a silent failure, and the reason
`mean_pairwise_identity` is reported beside every verdict. A block above ~95%
identity is untested, not negative, and is labelled `underpowered` rather than
`no_covariation`.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_ncrna.compute import require_slurm
from protein_ncrna.seqio import read_fasta

UNDERPOWERED_ID = 0.95          # fallback only, when R-scape power is absent
MIN_EXPECTED_PAIRS = 1.0        # R-scape expected-to-covary floor


def have(name: str) -> str | None:
    return shutil.which(name)


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        cols = next(fh).rstrip("\n").split("\t")
        return [dict(zip(cols, line.rstrip("\n").split("\t"))) for line in fh]


def mean_identity(seqs: list[str], cap: int = 60) -> float:
    """Mean pairwise identity over aligned columns, on a capped sample."""
    s = seqs[:cap]
    if len(s) < 2:
        return 1.0
    tot, n = 0.0, 0
    for i in range(len(s)):
        for j in range(i + 1, len(s)):
            a, b = s[i], s[j]
            cols = [(x, y) for x, y in zip(a, b) if x != "-" or y != "-"]
            if not cols:
                continue
            tot += sum(1 for x, y in cols if x == y) / len(cols)
            n += 1
    return tot / n if n else 1.0


def run_mafft(fa: Path, out: Path, threads: int) -> bool:
    exe = have("mafft")
    if not exe:
        return False
    proc = subprocess.run([exe, "--quiet", "--thread", str(threads), "--auto", str(fa)],
                          stdout=open(out, "w"), stderr=subprocess.PIPE)
    return proc.returncode == 0 and out.stat().st_size > 0


def to_stockholm(aln: Path, sto: Path, ss_cons: str = "") -> int:
    """Stockholm for R-scape, with the RNAalifold consensus as #=GC SS_cons.

    Two things R-scape needs that a naive dump does not give it:

      * a proposed structure. Without SS_cons, `R-scape -s` refuses with
        "Nucleotide alignment does not include a structure" -- the two-set test
        compares covariation inside the proposed helices against everything
        else, so it has nothing to work with.
      * short sequence names. Anchor ids run past 60 characters and carry '|'
        and ':'; Easel's Stockholm parser chokes on them. Names are replaced
        with seqNNNN and the mapping is written beside the file so every row
        remains traceable to its locus.
    """
    recs = [(h.split()[0], s) for h, s in read_fasta(str(aln))]
    if not recs:
        return 0
    names = [f"seq{i:05d}" for i in range(len(recs))]
    w = max(len(n) for n in names)
    with open(sto.with_suffix(".names.tsv"), "w") as fh:
        fh.write("stockholm_name\tanchor_id\n")
        for n, (h, _) in zip(names, recs):
            fh.write(f"{n}\t{h}\n")
    with open(sto, "w") as fh:
        fh.write("# STOCKHOLM 1.0\n\n")
        for n, (_, s) in zip(names, recs):
            fh.write(f"{n.ljust(w)}  {s.upper().replace('T', 'U')}\n")
        if ss_cons and len(ss_cons) == len(recs[0][1]):
            # R-scape wants WUSS; alifold emits plain brackets, which are valid.
            fh.write(f"{'#=GC SS_cons'.ljust(w)}  {ss_cons}\n")
        fh.write("//\n")
    return len(recs)


def run_alifold(aln: Path, outdir: Path) -> dict:
    exe = have("RNAalifold")
    if not exe:
        return {"status": "RNAalifold_missing"}
    proc = subprocess.run([exe, "--noPS", "--input-format=F", str(aln)],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          cwd=str(outdir))
    if proc.returncode != 0:
        return {"status": f"alifold_failed:{proc.stderr.decode(errors='replace')[-80:]}"}
    txt = proc.stdout.decode(errors="replace")
    (outdir / "alifold.out").write_text(txt)
    lines = [l for l in txt.splitlines() if l.strip()]
    struct, energy, cov = "", None, None
    for l in lines:
        m = re.match(r"^([.()<>{}\[\],_-]+)\s+\(\s*(-?[\d.]+)\s*=\s*(-?[\d.]+)\s*\+\s*(-?[\d.]+)\)", l)
        if m:
            struct, energy, cov = m.group(1), float(m.group(2)), float(m.group(4))
            break
    paired = struct.count("(")
    return {"status": "ok", "structure": struct, "consensus_mfe": energy,
            "covariance_term": cov, "n_base_pairs": paired,
            "struct_len": len(struct)}


def run_rscape(sto: Path, outdir: Path, evalue: float) -> dict:
    """R-scape two-set test, parsed from its own report rather than guessed.

    Three numbers matter and all come from R-scape itself:
      * `nbpairs` / `observed_covarying` -- how many proposed pairs covary;
      * `expected_covarying` -- how many it *could* have detected given the
        substitutions actually present in the alignment. This is the honest
        power measure. When it is near zero the alignment is too conserved for
        covariation to show up, and a null result says nothing about structure.
    """
    exe = have("R-scape")
    if not exe:
        return {"status": "rscape_missing"}
    d = outdir / "rscape"
    # Start from an empty directory. The power file is found by glob, so a
    # stale one left by an earlier run of the same block would be read as this
    # run's power analysis -- and power is what separates `underpowered` from
    # `tested_no_covariation`, so a stale read changes the verdict silently.
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    proc = subprocess.run([exe, "-s", "--outdir", str(d), "-E", str(evalue), str(sto)],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    txt = proc.stdout.decode(errors="replace") + proc.stderr.decode(errors="replace")
    (d / "rscape.log").write_text(txt)
    if proc.returncode != 0:
        last = txt.strip().splitlines()[-1][:80] if txt.strip() else "?"
        return {"status": f"rscape_failed:{last}"}

    out: dict = {"status": "ok"}
    m = re.search(r"# MSA \S+ nseq (\d+).*?alen (\d+).*?avgid ([\d.]+).*?nbpairs (\d+)",
                  txt)
    if m:
        out.update(rscape_nseq=int(m.group(1)), rscape_alen=int(m.group(2)),
                   rscape_avgid=float(m.group(3)), bpairs_tested=int(m.group(4)))
    # "[FP | TP True Found | Sen PPV F]" on the method line.
    m = re.search(r"^#\s*GTp\s+\S+\s+\[[^\]]*\]\s+\[\s*(\d+)\s*\|\s*(\d+)\s+"
                  r"(\d+)\s+(\d+)\s*\|\s*([\d.]+)\s+([\d.]+)", txt, re.M)
    if m:
        # [FP | TP True Found | Sen PPV F]. The number that matters is TP:
        # pairs that covary *and* lie in the proposed structure. "Found" counts
        # every significant pair including ones outside it, so reading Found
        # turns a lone false positive into apparent support -- which is exactly
        # what it did for an IS110 block scoring [1 | 0 22 1 | 0 0 0].
        out.update(false_positive_pairs=int(m.group(1)),
                   covarying_pairs=int(m.group(2)),
                   bpairs_tested=int(m.group(3)),
                   found_pairs=int(m.group(4)),
                   sensitivity=float(m.group(5)), ppv=float(m.group(6)))
    powers = sorted(d.glob("*.power"))
    if len(powers) > 1:
        # The directory is wiped before each run, so this means R-scape itself
        # emitted several. Refuse to pick one at random.
        out["parse_warning"] = (f"{len(powers)} .power files: "
                                f"{[p.name for p in powers]}")
    for pf in powers[:1]:
        ptxt = pf.read_text()
        m = re.search(r"# avg substitutions per BP\s+([\d.]+)", ptxt)
        if m:
            out["avg_subs_per_bp"] = float(m.group(1))
        m = re.search(r"# BPAIRS expected to covary\s+([\d.]+)", ptxt)
        if m:
            out["expected_covarying"] = float(m.group(1))
        m = re.search(r"# BPAIRS observed to covary\s+(\d+)", ptxt)
        if m:
            out["observed_covarying"] = int(m.group(1))

    out.setdefault("covarying_pairs", out.get("observed_covarying"))
    # The power file reports the same quantity independently; disagreement means
    # the report was misparsed, so surface it rather than trusting one side.
    obs = out.get("observed_covarying")
    if obs is not None and out.get("covarying_pairs") not in (None, obs):
        out["parse_warning"] = (f"GTp TP={out['covarying_pairs']} != "
                                f"power observed={obs}")
        out["covarying_pairs"] = min(out["covarying_pairs"], obs)
    return out


COLS = ["block_id", "family", "clade_id", "score", "n_seqs", "aln_len",
        "mean_pairwise_identity", "rscape_avgid", "consensus_mfe",
        "covariance_term", "n_base_pairs", "bpairs_tested", "covarying_pairs",
        "expected_covarying", "found_pairs", "false_positive_pairs",
        "avg_subs_per_bp", "sensitivity", "ppv",
        "verdict", "alifold_status", "rscape_status", "structure_file"]


# A verdict of `supported` requires all three: enough power for the test to
# have meant something, enough covarying pairs that one lucky pair cannot carry
# it, and most significant pairs landing inside the proposed structure.
SUPPORTED_MIN_EXPECTED = 3.0     # matches the triage `adequate` power band
SUPPORTED_MIN_TP = 3
SUPPORTED_MIN_PPV = 0.80


def verdict(ident: float, ali: dict, rs: dict) -> str:
    """Separate "no structure" from "could not have seen one", in four grades.

    Power comes from R-scape's own analysis -- the number of proposed pairs it
    expected to covary given the substitutions present -- not from a sequence
    identity rule of thumb. Below `MIN_EXPECTED_PAIRS` the test had no chance,
    so a null is reported as `underpowered`, which is the silent-failure mode
    this whole pipeline is built to avoid.

    The earlier rule promoted any block with a single covarying pair. That let a
    block with TP 1 against FP 2 and PPV 33% read as `supported`, which is a
    labelling error, not a marginal call: at expected-to-covary 0.1 a lone pair
    is what the null produces. `supported` now needs adequate power *and* three
    pairs *and* 80% PPV; anything with a real pair that misses those is
    `borderline`, which keeps the signal visible without overstating it.

    `tested_no_covariation` is the only genuinely negative verdict, and it is
    reserved for alignments that had the power to show covariation and did not.
    """
    if ali.get("status") != "ok":
        return "fold_failed"
    if rs.get("status") == "rscape_missing":
        return "structure_only_no_covariation_test"
    if rs.get("status", "").startswith("rscape_failed"):
        return "covariation_test_failed"
    cp = rs.get("covarying_pairs")
    if cp is None:
        return "covariation_test_failed"
    exp = rs.get("expected_covarying")
    ppv = rs.get("ppv")
    if (exp is not None and exp >= SUPPORTED_MIN_EXPECTED
            and cp >= SUPPORTED_MIN_TP
            and ppv is not None and ppv / 100.0 >= SUPPORTED_MIN_PPV):
        return "covariation_supported"
    if cp > 0:
        return "covariation_borderline"
    if exp is not None and exp >= SUPPORTED_MIN_EXPECTED:
        return "tested_no_covariation"
    if exp is not None:
        return "underpowered"
    # No power section at all: R-scape omits it when the alignment has no
    # substitutions, which is the extreme of underpowered, not a missing test.
    return "underpowered"


def reclassify(path: Path) -> int:
    """Rewrite the verdict column of an existing fold_summary.tsv in place.

    Re-running MAFFT and R-scape to change a label would be wasteful and would
    also make the old and new calls incomparable, since R-scape is stochastic in
    its null sampling. Every input to `verdict` is already a column here.
    """
    with open(path) as fh:
        cols = next(fh).rstrip("\n").split("\t")
        rows = [dict(zip(cols, l.rstrip("\n").split("\t"))) for l in fh]
    f = lambda r, k: (None if r.get(k, "") == "" else float(r[k]))
    changed = []
    for r in rows:
        rs = {"status": r.get("rscape_status", ""),
              "covarying_pairs": f(r, "covarying_pairs"),
              "expected_covarying": f(r, "expected_covarying"),
              "ppv": f(r, "ppv")}
        new = verdict(f(r, "mean_pairwise_identity") or 0.0,
                      {"status": r.get("alifold_status", "")}, rs)
        if new != r["verdict"]:
            changed.append((r["block_id"], r["verdict"], new))
        r["verdict"] = new
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    v = defaultdict(int)
    for r in rows:
        v[r["verdict"]] += 1
    j = path.parent / "fold_summary.json"
    summary = json.loads(j.read_text()) if j.exists() else {}
    summary["verdicts"] = dict(sorted(v.items()))
    summary["verdict_rule"] = {
        "supported_requires": {"expected_covarying": SUPPORTED_MIN_EXPECTED,
                               "TP": SUPPORTED_MIN_TP, "PPV": SUPPORTED_MIN_PPV},
        "note": "Relabelled from existing R-scape output; no realignment.",
    }
    j.write_text(json.dumps(summary, indent=2))
    for b, o, n in changed:
        print(f"  {b[-34:]}  {o} -> {n}")
    print(f"\n{len(changed)} of {len(rows)} verdicts changed")
    print(json.dumps(dict(sorted(v.items())), indent=2))
    print(f"-> {path}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", type=Path,
                    help="export_blocks.py output dir")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--reclassify", type=Path,
                    help="existing fold_summary.tsv; rewrite its verdict column "
                         "from the R-scape numbers already in it, no realignment")
    ap.add_argument("--which", default="padded", choices=["padded", "core"])
    ap.add_argument("--min-score", type=float, default=0.0)
    ap.add_argument("--max-seqs", type=int, default=200,
                    help="cap sequences per alignment; a clade of 500 near-"
                         "identical members adds cost, not covariation power")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--evalue", type=float, default=0.05)
    args = ap.parse_args()
    require_slurm("step 6 structure test")
    if args.reclassify:
        return reclassify(args.reclassify)
    if not (args.fasta and args.out):
        ap.error("--fasta and --out are required unless --reclassify is given")
    args.out.mkdir(parents=True, exist_ok=True)

    for tool in ("mafft", "RNAalifold", "R-scape"):
        print(f"  {tool:<12} {have(tool) or 'NOT FOUND'}")

    pad = read_tsv(args.fasta / "block_padding.tsv")
    # The candidate set is fixed by export_blocks.py (which holds the
    # benchmark_only flag); this only applies a further score floor.
    todo = [r for r in pad if float(r["score"]) >= args.min_score]
    print(f"\n{len(todo)} blocks to fold ({args.which})")

    rows = []
    for i, r in enumerate(todo, 1):
        bid = r["block_id"]
        d = args.out / bid
        d.mkdir(parents=True, exist_ok=True)
        src = Path(r[f"{args.which}_fasta"])
        recs = [(h, s) for h, s in read_fasta(str(src))][:args.max_seqs]
        sub = d / "input.fa"
        with open(sub, "w") as fh:
            for h, s in recs:
                fh.write(f">{h}\n{s}\n")

        row = {c: "" for c in COLS}
        row.update(block_id=bid, family=r["family"], clade_id=r["clade_id"],
                   score=r["score"], n_seqs=len(recs))
        aln = d / "aln.fa"
        if len(recs) < 3:
            row["verdict"] = "too_few_sequences"
            rows.append(row)
            continue
        if not run_mafft(sub, aln, args.threads):
            row["verdict"] = "align_failed"
            rows.append(row)
            continue
        aseqs = [s for _, s in read_fasta(str(aln))]
        ident = mean_identity(aseqs)
        row["aln_len"] = len(aseqs[0]) if aseqs else 0
        row["mean_pairwise_identity"] = round(ident, 4)

        ali = run_alifold(aln, d)
        if ali.get("structure"):
            (d / "structure.txt").write_text(
                f"# {bid}\n# consensus MFE {ali.get('consensus_mfe')} "
                f"kcal/mol, covariance term {ali.get('covariance_term')}\n"
                f"{ali['structure']}\n")
        sto = d / "aln.sto"
        to_stockholm(aln, sto, ali.get("structure", ""))
        rs = run_rscape(sto, d, args.evalue)

        row.update(consensus_mfe=ali.get("consensus_mfe", ""),
                   covariance_term=ali.get("covariance_term", ""),
                   n_base_pairs=ali.get("n_base_pairs", ""),
                   rscape_avgid=rs.get("rscape_avgid", ""),
                   covarying_pairs=rs.get("covarying_pairs", ""),
                   bpairs_tested=rs.get("bpairs_tested", ""),
                   expected_covarying=rs.get("expected_covarying", ""),
                   found_pairs=rs.get("found_pairs", ""),
                   false_positive_pairs=rs.get("false_positive_pairs", ""),
                   sensitivity=rs.get("sensitivity", ""),
                   ppv=rs.get("ppv", ""),
                   avg_subs_per_bp=rs.get("avg_subs_per_bp", ""),
                   alifold_status=ali.get("status", ""),
                   rscape_status=rs.get("status", ""),
                   verdict=verdict(ident, ali, rs),
                   structure_file=str(d / "structure.txt"))
        rows.append(row)
        print(f"  [{i}/{len(todo)}] {bid[:54]} id={ident:.2f} "
              f"{row['verdict']}", flush=True)

    with open(args.out / "fold_summary.tsv", "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in sorted(rows, key=lambda r: -float(r["score"] or 0)):
            fh.write("\t".join(str(r.get(c, "")) for c in COLS) + "\n")

    v = defaultdict(int)
    for r in rows:
        v[r["verdict"]] += 1
    ids = [float(r["mean_pairwise_identity"]) for r in rows
           if r["mean_pairwise_identity"] != ""]
    summary = {"n_blocks": len(rows), "verdicts": dict(v),
               "median_mean_pairwise_identity": (round(statistics.median(ids), 4)
                                                 if ids else None),
               "tools": {t: have(t) for t in ("mafft", "RNAalifold", "R-scape")},
               "note": ("A block above %.0f%% mean identity is reported as "
                        "underpowered, not negative: there were no substitutions "
                        "for covariation to appear in." % (UNDERPOWERED_ID * 100))}
    (args.out / "fold_summary.json").write_text(json.dumps(summary, indent=2))
    print("\n" + json.dumps(summary, indent=2))
    print(f"-> {args.out}/fold_summary.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
