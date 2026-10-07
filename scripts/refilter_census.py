#!/usr/bin/env python
"""Re-apply the current anchor config to an existing census, without re-running it.

The 2026-10-06 census ran before the config was corrected: TnpB was labelled
Cas12, and per-family length floors did not exist. Re-running would repeat ~8
core-hours of prodigal for no new information, since the expensive part (ORF
calling and HMM search) does not change — only the labelling rules do.

This reads `configs/blind/anchors.json` as the single source of truth and rewrites a
census directory with current labels and floors. Logic is not duplicated: model
to family comes from the config, and the length floor is the same
`min_aa_len` lookup `find_anchors` applies.

    refilter_census.py --census <census_dir> --out <corrected_dir>
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_ncrna.compute import require_slurm
from protein_ncrna.anchors import load_config
from protein_ncrna.seqio import read_fasta


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--config", type=Path,
                    default=Path(__file__).resolve().parents[1] / "configs/blind/anchors.json")
    args = ap.parse_args()
    require_slurm("census refilter")

    cfg = load_config(args.config)
    model_family = {m["name"]: fam for fam, spec in cfg["families"].items()
                    for m in spec["models"]}
    global_min = cfg["thresholds"]["min_aa_len"]
    fam_min = {fam: spec.get("min_aa_len", global_min)
               for fam, spec in cfg["families"].items()}

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "proteins").mkdir(exist_ok=True)

    # Anchor sequences, keyed by anchor_id, from the census per-family FASTAs.
    seq_of: dict[str, str] = {}
    for faa in sorted((args.census / "proteins").glob("*.faa")):
        for h, s in read_fasta(str(faa)):
            seq_of[h.split()[0]] = s
    print(f"loaded {len(seq_of)} anchor sequences")

    at = args.census / "anchors.tsv"
    with open(at) as fh:
        header = fh.readline().rstrip("\n")
        cols = {c: i for i, c in enumerate(header.split("\t"))}

    dropped = Counter()
    relabelled = Counter()
    # 537 assemblies contain every contig twice, so the same locus appears as
    # two identical anchor rows. One locus is one anchor.
    seen_ids: set[str] = set()
    fam_stats: dict[str, dict] = defaultdict(
        lambda: {"anchors": 0, "genomes": set(), "contigs": set(), "species": set(),
                 "complete_windows": 0, "partial_anchors": 0, "ambiguous": 0})
    faa_out: dict[str, object] = {}
    kept = total = 0

    with open(at) as fh, open(args.out / "anchors.tsv", "w") as out:
        next(fh)
        out.write(header + "\n")
        for line in fh:
            total += 1
            f = line.rstrip("\n").split("\t")
            old_fam, model = f[cols["family"]], f[cols["model"]]
            aid = f[cols["anchor_id"]]
            if aid in seen_ids:
                dropped["duplicate anchor_id (doubled contig in assembly)"] += 1
                continue
            seen_ids.add(aid)

            # The winning model decides the family under the current config.
            new_fam = model_family.get(model)
            if new_fam is None:
                dropped[f"model no longer configured: {model}"] += 1
                continue
            seq = seq_of.get(aid)
            if seq is None:
                dropped["no sequence found"] += 1
                continue
            if len(seq) < fam_min[new_fam]:
                dropped[f"{new_fam}: shorter than min_aa_len {fam_min[new_fam]}"] += 1
                continue
            if new_fam != old_fam:
                relabelled[f"{old_fam} -> {new_fam}"] += 1
                f[cols["family"]] = new_fam
                f[cols["all_families"]] = ",".join(
                    new_fam if x == old_fam else x
                    for x in f[cols["all_families"]].split(","))

            out.write("\t".join(f) + "\n")
            kept += 1

            s = fam_stats[new_fam]
            s["anchors"] += 1
            s["genomes"].add(f[cols["genome"]])
            s["contigs"].add(f"{f[cols['genome']]}|{f[cols['contig']]}")
            s["species"].add(f[cols["species"]])
            s["complete_windows"] += f[cols["window_complete"]] == "True"
            s["partial_anchors"] += f[cols["anchor_partial"]] == "True"
            s["ambiguous"] += f[cols["ambiguous"]] == "True"

            if new_fam not in faa_out:
                faa_out[new_fam] = open(args.out / "proteins" / f"{new_fam}.faa", "w")
            faa_out[new_fam].write(f">{aid} {new_fam} {model} "
                                   f"species={f[cols['species']]}\n{seq}\n")
    for fh_ in faa_out.values():
        fh_.close()

    with open(args.out / "summary_by_anchor_family.tsv", "w") as fh:
        fh.write("family\tanchors\tgenomes\tcontigs\tspecies\tcomplete_window_frac\t"
                 "partial_anchor_frac\tambiguous_frac\tanchors_per_genome\n")
        for fam, s in sorted(fam_stats.items(), key=lambda kv: -kv[1]["anchors"]):
            n = s["anchors"]
            fh.write(f"{fam}\t{n}\t{len(s['genomes'])}\t{len(s['contigs'])}\t"
                     f"{len(s['species'])}\t{s['complete_windows']/n:.4f}\t"
                     f"{s['partial_anchors']/n:.4f}\t{s['ambiguous']/n:.4f}\t"
                     f"{n/max(len(s['genomes']),1):.3f}\n")

    (args.out / "refilter_report.json").write_text(json.dumps({
        "source_census": str(args.census),
        "config": str(args.config),
        "anchors_in": total, "anchors_kept": kept, "anchors_dropped": total - kept,
        "dropped_reasons": dict(dropped.most_common()),
        "relabelled": dict(relabelled.most_common()),
        "per_family": {f: {"anchors": s["anchors"], "genomes": len(s["genomes"]),
                           "species": len(s["species"])} for f, s in fam_stats.items()},
    }, indent=2))

    print(f"kept {kept}/{total} anchors ({total - kept} dropped)")
    for k, v in relabelled.most_common():
        print(f"  relabelled {v:>7}  {k}")
    for k, v in dropped.most_common():
        print(f"  dropped    {v:>7}  {k}")
    print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
