#!/usr/bin/env python
"""One annotated GenBank record per case, one subfolder per bag.

    export_loci_gbk.py --windows <dir ...> --blocks <blocks.tsv> \\
                       --members <block_members.tsv> --census <anchors.tsv> \\
                       --out <dir> [--bags <list>] [--ranked <tsv>] \\
                       [--pfam <tsv>] [--pfam-sensitive <tsv>] \\
                       [--orig <insertions_bags.tsv>]

A case is one genomic locus: one anchor protein hit, with the window carved
around it. Every case in this corpus is a tblastn extension -- the anchors come
from expand_by_homology against all real genomes, not from the original
insertion-pair extraction -- so "the bag and its BLAST-extended cases" is the
whole window set, and each case gets a record whether or not a block was called
in it. A case with no ncRNA feature is a result too: it says this locus carries
the protein and not the RNA.

Coordinates. The stored window is already in ANCHOR orientation: the anchor ORF
runs forward at [anchor_offset, anchor_offset + anchor_len) regardless of its
strand on the contig. Block member coordinates are in the anchor frame,
x = window_pos - anchor_offset, so a feature lands at anchor_offset + x. Every
record is written in that frame, which is the frame the blocks were scored in.

Genomic coordinates are also emitted, but only when they can be proven. The
anchor's true contig span is embedded in its own id, so the window->contig
mapping is applied to the anchor first and checked against it; if it does not
round-trip, this case's genomic coordinates are left blank rather than guessed.
Window geometry is the one thing here that has bitten this pipeline before.

No ORFs beyond the anchor. Step 4's coding mask was not stored, so predicting
genes again would put a second, slightly different ORF set next to blocks whose
edges were defined by the first one. An ORF in this file that did not bound the
block it abuts would be read as if it had.

Post-hoc: this renders results that already exist. It is not in BLIND_CHAIN and
must never feed anything that selects or scores.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# anchor_id looks like GCA_003734765.1|RKXS01000016.1:84543-85766:-
# and the trailing span is the anchor's true 1-based inclusive contig span.
ANCHOR_SPAN = re.compile(r":(\d+)-(\d+):([+-])$")


def case_name(anchor_id: str) -> str:
    """A filename that can be mapped back to the anchor it came from."""
    s = anchor_id.replace("|", "__").replace(":", "_")
    s = s.replace("+", "plus").replace("-", "-")
    return re.sub(r"[^A-Za-z0-9._-]", "_", s)


def genomic(p0: int, p1: int, win_start: int, win_end: int, strand: int):
    """Window half-open [p0,p1) -> contig 1-based inclusive (start, end).

    win_start and win_end are themselves 1-based inclusive: their span is
    win_end - win_start + 1, which equals the stored sequence length. Treating
    win_start as 0-based shifts every plus-strand feature one base left, which
    is invisible in a length check and wrong in every record.
    """
    if strand == -1:
        return win_end - p1 + 1, win_end - p0
    return win_start + p0, win_start + p1 - 1


def f(x, d=""):
    return x if x not in (None, "") else d


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", required=True, type=Path, nargs="+")
    ap.add_argument("--blocks", required=True, type=Path)
    ap.add_argument("--members", required=True, type=Path)
    ap.add_argument("--census", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--bags", type=Path,
                    help="bag list; default is every bag with a real block")
    ap.add_argument("--ranked", type=Path)
    ap.add_argument("--candidates", type=Path)
    ap.add_argument("--pfam", type=Path)
    ap.add_argument("--pfam-sensitive", type=Path)
    ap.add_argument("--orig", type=Path,
                    help="bags_v4 insertions_bags.tsv, for the bag's original "
                         "insertion-site statistics")
    ap.add_argument("--date", default="09-OCT-2026",
                    help="GenBank LOCUS date; fixed so reruns are byte-stable")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    args = ap.parse_args()

    blocks = [r for r in csv.DictReader(open(args.blocks), delimiter="\t")
              if r["decoy"] == "False"]
    by_bag_blocks = defaultdict(list)
    for r in blocks:
        by_bag_blocks[r["clade"]].append(r)

    bags = (args.bags.read_text().split() if args.bags
            else sorted(by_bag_blocks))
    bags = [b for i, b in enumerate(sorted(bags))
            if i % args.nshards == args.shard]
    bagset = set(bags)

    # members: bag -> anchor_id -> [block rows]
    mem = defaultdict(lambda: defaultdict(list))
    for r in csv.DictReader(open(args.members), delimiter="\t"):
        if r["decoy"] == "False" and r["clade"] in bagset:
            mem[r["clade"]][r["anchor_id"]].append(r)

    # Side tables, all optional; a missing one must leave a blank field rather
    # than drop the bag.
    def load(p, key):
        if not p or not p.exists():
            return {}
        return {r[key]: r for r in csv.DictReader(open(p), delimiter="\t")}

    ranked = load(args.ranked, "bag")
    cands = load(args.candidates, "bag")
    pfam = load(args.pfam, "bag")
    psens = load(args.pfam_sensitive, "bag")
    orig = load(args.orig, "bag_id")

    # The census is far too large to hold entire; keep only the selected bags.
    cen = {}
    if args.census and args.census.exists():
        for r in csv.DictReader(open(args.census), delimiter="\t"):
            if r.get("family") in bagset:
                cen[r["anchor_id"]] = r

    src_of = {}
    for root in args.windows:
        if root.is_dir():
            for d in sorted(root.iterdir()):
                if d.is_dir() and (d / f"{d.name}.jsonl").exists():
                    src_of.setdefault(d.name, d)

    args.out.mkdir(parents=True, exist_ok=True)
    tot_cases = tot_rna = tot_nogeo = 0
    for bag in bags:
        d = src_of.get(bag)
        if d is None:
            print(f"  {bag}: no windows, skipped")
            continue
        bd = args.out / bag
        (bd / "gbk").mkdir(parents=True, exist_ok=True)

        bb = sorted(by_bag_blocks.get(bag, []),
                    key=lambda r: -float(r["composite_score"]))
        bcols = ["block_id", "composite_score", "recurrence", "member_fraction",
                 "n_members", "n_clade_windows", "median_start", "median_len",
                 "dist_to_anchor", "side", "coord_tightness",
                 "strand_consistency", "length_consistency", "median_pident"]
        with open(bd / "blocks.tsv", "w") as fh:
            fh.write("\t".join(bcols) + "\n")
            for r in bb:
                fh.write("\t".join(str(f(r.get(c))) for c in bcols) + "\n")

        cases, n_rna = [], 0
        for line in open(d / f"{bag}.jsonl"):
            w = json.loads(line)
            aid = w["anchor_id"]
            seq = w["seq"]
            ao, al = w["anchor_offset"], w["anchor_len"]
            ws, we, st = w["win_start"], w["win_end"], int(w["strand"])
            c = cen.get(aid, {})

            # Prove the window->contig mapping on the anchor before trusting it
            # for anything else.
            m = ANCHOR_SPAN.search(aid)
            gs, ge = genomic(ao, ao + al, ws, we, st)
            geo_ok = bool(m) and (int(m.group(1)), int(m.group(2))) == (gs, ge)
            if not geo_ok:
                tot_nogeo += 1

            # The LOCUS name must be unique per case and short; the full
            # anchor_id is too long for it and is carried in source/note and in
            # the filename instead, so nothing depends on a truncated string.
            rec = SeqRecord(
                Seq(seq), id=f"{bag}_{len(cases):05d}",
                name=f"{bag}_{len(cases):05d}"[:16],
                description=(f"{bag} anchor-tethered locus, "
                             f"{f(w.get('species'), 'unknown')}"),
                annotations={
                    "molecule_type": "DNA", "topology": "linear",
                    "data_file_division": "BCT",
                    "date": args.date,
                    "organism": f(c.get("organism"),
                                  f(w.get("species"), "unknown")),
                    "source": f(c.get("organism"),
                                f(w.get("species"), "unknown")),
                    "comment": (
                        "Blind protein-anchored ncRNA discovery. Coordinates "
                        "are in ANCHOR orientation: the anchor ORF runs forward "
                        "and all features are placed at anchor_offset + x, "
                        "where x is the anchor-frame coordinate the blocks were "
                        "scored in. This case is a tblastn homology extension "
                        "of the bag representative. Only the anchor CDS and "
                        "called ncRNA blocks are annotated; step 4's coding "
                        "mask was not retained, so no other ORFs are drawn. "
                        + ("Genomic coordinates verified against the anchor's "
                           "own contig span." if geo_ok else
                           "Genomic coordinates OMITTED: the window-to-contig "
                           "mapping did not round-trip on this anchor.")),
                })

            astr = {"1": "+", "-1": "-", "+": "+", "-": "-"}.get(
                str(c.get("anchor_strand", "")), "?")
            note = [f"anchor_id {aid}",
                    f"window {ws}-{we} on {f(w.get('contig'))}",
                    f"window_strand {'+' if st == 1 else '-'}",
                    f"anchor_offset {ao}",
                    "origin tblastn expand_by_homology"]
            rec.features.append(SeqFeature(
                SimpleLocation(0, len(seq), strand=1), type="source",
                qualifiers={
                    "organism": [f(c.get("organism"),
                                   f(w.get("species"), "unknown"))],
                    "mol_type": ["genomic DNA"],
                    "db_xref": [f"BioSample:{f(c.get('genome'))}"]
                    if c.get("genome") else [],
                    "note": ["; ".join(note)],
                    "chromosome": [f(w.get("contig"))],
                }))

            aq = {
                "locus_tag": [bag], "gene": [bag],
                "product": ["anchor protein (bag representative homolog)"],
                "note": [f"bag {bag}; tblastn bitscore "
                         f"{f(c.get('score'), 'NA')}; evalue "
                         f"{f(c.get('evalue'), 'NA')}; contig strand {astr}"],
            }
            if geo_ok:
                aq["note"][0] += f"; contig {gs}..{ge}"
            rec.features.append(SeqFeature(
                SimpleLocation(ao, ao + al, strand=1), type="CDS",
                qualifiers=aq))

            hits = mem[bag].get(aid, [])
            for h in sorted(hits, key=lambda r: int(r["start"])):
                p0 = ao + int(r_start := int(h["start"]))
                p1 = ao + int(h["end"])
                p0c, p1c = max(0, p0), min(len(seq), p1)
                if p1c <= p0c:
                    continue
                blk = next((r for r in bb if r["block_id"] == h["block_id"]),
                           {})
                q = {
                    "ncRNA_class": ["other"],
                    "product": [f"candidate ncRNA block {h['block_id']}"],
                    "note": [
                        f"block {h['block_id']}; anchor_frame "
                        f"{r_start}..{h['end']}; side {f(blk.get('side'))}; "
                        f"dist_to_anchor {f(blk.get('dist_to_anchor'))}; "
                        f"composite_score {f(blk.get('composite_score'))}; "
                        f"recurrence {f(blk.get('recurrence'))}; "
                        f"member_fraction {f(blk.get('member_fraction'))}; "
                        f"pident_this_case {h['pident']}; "
                        f"blind call, no RNA family model used"],
                }
                if geo_ok:
                    bgs, bge = genomic(p0c, p1c, ws, we, st)
                    q["note"][0] += f"; contig {bgs}..{bge}"
                rec.features.append(SeqFeature(
                    SimpleLocation(p0c, p1c, strand=1), type="ncRNA",
                    qualifiers=q))
            if hits:
                n_rna += 1

            with open(bd / "gbk" / f"{case_name(aid)}.gbk", "w") as fh:
                SeqIO.write(rec, fh, "genbank")

            cases.append({
                "case_file": f"gbk/{case_name(aid)}.gbk",
                "anchor_id": aid, "genome": f(c.get("genome")),
                "species": f(w.get("species")),
                "organism": f(c.get("organism")),
                "contig": f(w.get("contig")),
                "win_start": ws, "win_end": we,
                "window_strand": "+" if st == 1 else "-",
                "anchor_offset": ao, "anchor_len": al,
                "window_len": len(seq),
                "window_complete": w.get("window_complete", ""),
                "tblastn_bitscore": f(c.get("score")),
                "tblastn_evalue": f(c.get("evalue")),
                "genomic_coords_verified": geo_ok,
                "n_ncrna_features": len(hits),
                "block_ids": ";".join(sorted({h["block_id"] for h in hits})),
            })

        ccols = list(cases[0]) if cases else []
        with open(bd / "cases.tsv", "w") as fh:
            fh.write("\t".join(ccols) + "\n")
            for r in cases:
                fh.write("\t".join(str(r[c]) for c in ccols) + "\n")

        rk, cd = ranked.get(bag, {}), cands.get(bag, {})
        og = orig.get(bag, {})
        info = {
            "bag": bag,
            "n_cases": len(cases),
            "n_cases_with_an_ncrna_block": n_rna,
            "fraction_of_cases_with_a_block": (
                round(n_rna / len(cases), 4) if cases else 0.0),
            "n_blocks": len(bb),
            "max_composite_score": (
                float(bb[0]["composite_score"]) if bb else None),
            "case_origin": "tblastn expand_by_homology vs all real genomes",
            "coordinate_frame": "anchor orientation; x = window_pos - "
                                "anchor_offset",
            "n_cases_without_verified_genomic_coords": sum(
                1 for c in cases if not c["genomic_coords_verified"]),
            "lane": f(cd.get("lane")),
            "n_distinct_anchors": f(cd.get("n_distinct_anchors")),
            "pfam_curated_domains": f(pfam.get(bag, {}).get("all_domains")),
            "pfam_sensitive_status": f(
                psens.get(bag, {}).get("sensitive_status")),
            "pfam_sensitive_leads": f(
                psens.get(bag, {}).get("sensitive_leads")),
            "hard_null": {
                "rank": f(rk.get("rank")), "tier": f(rk.get("tier")),
                "evidence_score": f(rk.get("evidence_score")),
                "potential_score": f(rk.get("potential_score")),
                "information_gap": f(rk.get("information_gap")),
                "n_decoy_clades": f(rk.get("n_decoy_clades")),
                "n_decoy_blocks": f(rk.get("n_decoy_blocks")),
                "margin": f(rk.get("margin")),
                "challenge_basis": f(rk.get("challenge_basis")),
                "tier_reason": f(rk.get("tier_reason")),
            } if rk else "not among the hard-null targets",
            "original_bags_v4": {
                "n_sites": f(og.get("n_sites")),
                "n_species": f(og.get("n_species")),
                "n_phyla": f(og.get("n_phyla")),
                "median_nc_total_len": f(og.get("median_nc_total_len")),
                "species": f(og.get("species")),
            } if og else "not in insertions_bags.tsv",
        }
        (bd / "bag_info.json").write_text(json.dumps(info, indent=2))

        rl = [f"{bag}", "=" * len(bag), "",
              f"{len(cases)} cases, {n_rna} carrying a called ncRNA block "
              f"({info['fraction_of_cases_with_a_block']:.1%}).",
              f"{len(bb)} blocks; best composite score "
              f"{info['max_composite_score']}.", "",
              "Every case is a tblastn homology extension of the bag",
              "representative, found against all real genomes. A case with no",
              "ncRNA feature carries the protein and not the RNA, which is a",
              "result and not a gap.", "",
              "GenBank coordinates are in anchor orientation: the anchor ORF",
              "runs forward and blocks sit at anchor_offset + x. Genomic",
              "coordinates appear only where the window-to-contig mapping was",
              "verified against the anchor's own span.", "",
              "Only the anchor CDS and called blocks are annotated. Step 4's",
              "coding mask was not retained, so no other ORFs are drawn rather",
              "than drawing a second set that did not bound these blocks.", "",
              f"Pfam (curated): {info['pfam_curated_domains'] or 'none'}",
              f"Pfam (sub-threshold): "
              f"{info['pfam_sensitive_leads'] or info['pfam_sensitive_status'] or 'n/a'}",
              ""]
        if rk:
            rl += [f"Hard null: tier {rk.get('tier')}, rank {rk.get('rank')} "
                   f"of 18. {rk.get('tier_reason')}.",
                   f"  evidence {rk.get('evidence_score')} / potential "
                   f"{rk.get('potential_score')} / gap "
                   f"{rk.get('information_gap')}",
                   f"  margin {rk.get('margin')} vs "
                   f"{rk.get('challenge_basis')}, from "
                   f"{rk.get('n_decoy_clades')} decoy clades producing "
                   f"{rk.get('n_decoy_blocks')} blocks", ""]
        (bd / "README.txt").write_text("\n".join(rl))

        tot_cases += len(cases)
        tot_rna += n_rna
        print(f"  {bag}: {len(cases):5d} cases, {n_rna:5d} with a block, "
              f"{len(bb)} blocks")

    print(f"\n{len(bags)} bags, {tot_cases} cases, {tot_rna} with an ncRNA "
          f"feature, {tot_nogeo} without verified genomic coordinates")
    print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
