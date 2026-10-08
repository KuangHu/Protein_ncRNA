#!/usr/bin/env python
"""Convert an MGE insertion bag set into the window format step 4 consumes.

    bags_to_windows.py --bags <bags_v4> --cds <cds_v3> --inserts <insertions_inserts.fna> \\
                       --out <windows_dir> [--max-per-clade 500]

The insertion catalogue from `fna_based_mgefinder_project` is structurally the
same object this pipeline already works on. A genomic window is an anchor ORF
with flanking DNA in the anchor frame; an insertion is an anchor ORF (the
insert's dominant CDS) with the rest of the element around it. So the insert
*is* the window, with a boundary set by the element itself rather than by an
arbitrary +/-5 kb.

What that buys over the genomic corpus is the thing the genomic run ran out of:
distinct sequences. A 6,964-site bag carries 6,763 distinct non-coding
sequences, where 34 of 75 genomic clades had fewer than 20 distinct proteins and
were underpowered for covariation before folding began.

What it costs is window length. Inserts run ~0.8-1.3 kb against 11 kb windows,
so `coord_tightness` and `recurrence` are computed over a much shorter span and
the score floors derived from the genomic decoys do NOT transfer. The null is
rebuilt on this corpus and the floors must be re-derived from it.

Blind, on the same terms as the rest of the discovery chain. The clade is the
`cds_cluster_id` -- an MMseqs cluster of ab-initio Prodigal ORFs, with no HMM,
no IS library and no functional claim attached. The non-coding part is never
examined here; it is carried as sequence for step 4 to find blocks in.

Orientation: each insert is emitted on its dominant ORF's strand, so `x = 0` is
the anchor's first base and negative is upstream, exactly as in the genomic
windows. That matters because the gff gives half the dominant ORFs on the minus
strand, and leaving them would place the same element's non-coding region at
two mirrored coordinates and destroy the recurrence signal it is scored on.

Outputs the layout extract_windows.py produces: `<clade>.fasta`, `<clade>.jsonl`,
`proteins/<clade>.faa`, `cap_manifest.tsv`, `windows_summary.json`.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_ncrna.compute import require_slurm
from protein_ncrna.diversity import stratified_sample
from protein_ncrna.seqio import read_fasta

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

# Clade names are `<band>__<bag_id>`, because permutation_null.py reads the
# family as everything before the first `__` and draws each decoy from the
# *other* clades of that family. Naming the clades `CDS00002` would make every
# bag its own family with a pool of one clade, and the null would be empty while
# still reporting an FDR of 0 -- the same degenerate null that invalidated Cas3
# in the genomic run.
#
# The band is insert length, not an arbitrary grouping. The genomic null needed
# no length matching because every window was 11 kb; inserts run 500-5300 bp, and
# drawing a decoy for a 700 bp element out of 3 kb elements would make the decoy
# incoherent for reasons that have nothing to do with the relatedness the
# permutation is meant to break. An easier null inflates the real-vs-decoy
# separation, which is the dangerous direction to be wrong in.
BANDS = [(700, "L0lt700"), (1000, "L1_700_1000"), (1400, "L2_1000_1400"),
         (2000, "L3_1400_2000"), (10 ** 9, "L4ge2000")]


def band_of(median_len: float) -> str:
    for hi, name in BANDS:
        if median_len < hi:
            return name
    return BANDS[-1][1]


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def load_gff_dominant(gff: Path, want: dict[str, str]) -> dict[str, tuple]:
    """db_id -> (start0, end0, strand) for its dominant ORF.

    Prodigal writes one record per ORF with the insert as the sequence id and
    `ID=<n>_<m>` in the attributes; the dominant ORF id from insert_cds.tsv is
    `<db_id>_<m>`, so the ORF is matched on that suffix rather than on position.
    """
    out: dict[str, tuple] = {}
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        db = f[0]
        tgt = want.get(db)
        if tgt is None or db in out:
            continue
        # ID=1_2 -> this is ORF 2 of the insert; dominant_orf is "<db_id>_2".
        attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        n = attr.get("ID", "_").split("_")[-1]
        if f"{db}_{n}" != tgt:
            continue
        out[db] = (int(f[3]) - 1, int(f[4]), f[6])
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bags", required=True, type=Path)
    ap.add_argument("--cds", required=True, type=Path)
    ap.add_argument("--inserts", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--max-per-clade", type=int, default=500)
    ap.add_argument("--min-sites", type=int, default=8,
                    help="skip bags with fewer sites than this")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    require_slurm("bag window conversion")
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "proteins").mkdir(exist_ok=True)

    # --- which sites are in the bag set ----------------------------------
    bag_of: dict[str, str] = {}
    species: dict[str, str] = {}
    for line in open(args.bags / "insertions_sites.jsonl"):
        d = json.loads(line)
        db = d["site_id"][:-5] if d["site_id"].endswith(".site") else d["site_id"]
        bag_of[db] = d["bag_id"]
        species[db] = d.get("species", "")
    print(f"{len(bag_of)} sites in {len(set(bag_of.values()))} bags")

    # --- dominant ORF per insert -----------------------------------------
    dom: dict[str, str] = {}
    meta: dict[str, dict] = {}
    with open(args.cds / "insert_cds.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["db_id"] in bag_of:
                dom[r["db_id"]] = r["dominant_orf"]
                meta[r["db_id"]] = r
    print(f"{len(dom)} of those have a dominant ORF in insert_cds.tsv")

    orf = load_gff_dominant(args.cds / "orfs.gff", dom)
    print(f"{len(orf)} dominant ORFs located in orfs.gff")

    prot = {h.split()[0]: s for h, s in read_fasta(str(args.cds / "dominant.faa"))}
    print(f"{len(prot)} dominant proteins loaded")

    # --- insert sequences --------------------------------------------------
    seqs = {h.split()[0]: s for h, s in read_fasta(str(args.inserts))}
    print(f"{len(seqs)} insert sequences loaded")

    # --- group by bag, emit in the anchor frame ---------------------------
    groups: dict[str, list[str]] = defaultdict(list)
    for db in sorted(dom):
        if db in orf and db in seqs:
            groups[bag_of[db]].append(db)
    groups = {b: m for b, m in groups.items() if len(m) >= args.min_sites}
    print(f"{len(groups)} bags with >= {args.min_sites} usable sites, "
          f"{sum(len(m) for m in groups.values())} inserts")

    manifest = open(args.out / "cap_manifest.tsv", "w")
    manifest.write("clade_id\tanchor_id\tselected\treason\n")
    n_win = n_minus = 0
    summary = []
    for bag in sorted(groups):
        mem = sorted(groups[bag])
        pseq = {m: prot[dom[m]] for m in mem if dom[m] in prot}
        if len(mem) > args.max_per_clade:
            keep = set(stratified_sample(mem, pseq, args.max_per_clade,
                                         seed=args.seed))
            reason = "stratified_kmer"
        else:
            keep, reason = set(mem), "under_cap"
        med_len = statistics.median([len(seqs[m]) for m in keep])
        clade = f"{band_of(med_len)}__{bag}"
        for m in mem:
            manifest.write(f"{clade}\t{m}\t{m in keep}\t{reason}\n")

        fa = open(args.out / f"{clade}.fasta", "w")
        jl = open(args.out / f"{clade}.jsonl", "w")
        faa = open(args.out / "proteins" / f"{clade}.faa", "w")
        kept = 0
        for m in sorted(keep):
            s, e, strand = orf[m]
            seq = seqs[m]
            if strand == "-":
                n_minus += 1
                seq = rc(seq)
                s, e = len(seq) - e, len(seq) - s
            rec = {
                "anchor_id": m, "contig": m,
                "win_start": 0, "win_end": len(seq),
                "strand": 1, "anchor_offset": s, "anchor_len": e - s,
                "family": band_of(med_len), "species": species.get(m, ""),
                "window_complete": True, "seq": seq.upper(),
            }
            jl.write(json.dumps(rec) + "\n")
            # discover_blocks.py parses the anchor frame out of `key=value`
            # pairs in the FASTA header, NOT from the jsonl. A bare `>id` leaves
            # anchor_offset and anchor_len at 0, which silently puts every
            # coordinate back in raw window space (no negative x, so every block
            # reads as "downstream") and masks the anchor ORF as zero-length, so
            # it is never treated as coding.
            fa.write(f">{m} family={rec['family']} species={rec['species']} "
                     f"anchor_offset={s} anchor_len={e - s}\n{seq.upper()}\n")
            if dom[m] in prot:
                faa.write(f">{m}\n{prot[dom[m]]}\n")
            kept += 1
            n_win += 1
        fa.close(), jl.close(), faa.close()
        summary.append({"clade_id": clade, "bag_id": bag, "n_sites": len(mem),
                        "n_windows": kept, "cap_reason": reason,
                        "median_window_len": int(med_len),
                        "n_species": len({species.get(m, "") for m in keep})})
        print(f"  {clade}: {len(mem)} sites -> {kept} windows "
              f"(median {med_len:.0f} bp, {reason})", flush=True)
    manifest.close()

    with open(args.out / "windows_summary.tsv", "w") as fh:
        cols = ["clade_id", "bag_id", "n_sites", "n_windows", "cap_reason",
                "median_window_len", "n_species"]
        fh.write("\t".join(cols) + "\n")
        for r in summary:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    (args.out / "windows_summary.json").write_text(json.dumps({
        "source_bags": str(args.bags), "clades": len(groups), "windows": n_win,
        "min_sites": args.min_sites, "max_per_clade": args.max_per_clade,
        "reverse_complemented": n_minus,
        "note": ("Score floors from the genomic run do not transfer: inserts are "
                 "~1 kb against 11 kb windows, so recurrence and coordinate "
                 "tightness are measured over a different span. Re-derive the "
                 "floors from this corpus's own permutation decoys."),
    }, indent=2))
    print(f"\n{n_win} windows in {len(groups)} bags "
          f"({n_minus} reverse-complemented to the anchor strand) -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
