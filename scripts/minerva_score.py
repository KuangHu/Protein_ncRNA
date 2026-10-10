#!/usr/bin/env python
"""Step 5 proposer: Minerva base-pairing structure over anchor windows.

    minerva_score.py --windows <dir of per-bag window dirs> --bags <list> \\
                     --out <dir> [--blocks <blocks.tsv>] [--max-cases 20] \\
                     [--model gbrixi/minerva-mlm-8k] [--threshold 0.6]

Why this exists. The block caller measures RECURRENCE: an RNA is called because
the same noncoding segment keeps appearing at the same offset from the anchor
across many lineages. That is strong evidence and it is also the method's
boundary -- 45 of 399 bags have fewer than the 5 distinct anchors a block needs,
and 175 more were scored and produced nothing. For those bags recurrence is not
weak evidence, it is undefined, and no amount of extra corpus fixes it.

Minerva asks a different question of a SINGLE locus: does this noncoding segment
look like it base-pairs with itself? That is orthogonal to recurrence, which is
exactly why it is worth running -- and exactly why it cannot stand in for the
null. A confident hairpin is not evidence an RNA is real or functional; plenty
of bacterial intergenic DNA folds. Minerva proposes candidates. It is never the
final evidence, and nothing it emits may re-enter the block caller, the decoys
or the floor.

On the blind principle: Minerva-MLM is a self-supervised genome language model
over prokaryotic genomes, not a curated RNA-family model. It carries no msr/msd,
tracrRNA, bridge-RNA or guide-RNA template, so using it does not import the
known-RNA priors the discovery stage excludes. Its outputs land in columns
prefixed `minerva_` so no reader can mistake them for blind evidence.

Genes are called by Minerva's own Pyrodigal pass rather than taken from our
GenBank files, which annotate the anchor CDS alone. Minerva needs every CDS to
build its mixed amino-acid/nucleotide representation; handing it one CDS would
feed the flanking genes in as raw nucleotides and misrepresent the locus.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def split_mixed_tokens(ts: str) -> list[str]:
    """Mixed-token string -> token list, matching ``mixed_token_length``.

    The tokenizer is character-level and adds no special tokens, so every
    character is one token except the strand markers ``<+>``/``<->``, which are
    one token each. call_structures needs this list to tell a nucleotide region
    (lowercase) from a coding one (amino acids), and getting it wrong silently
    shifts every region boundary.
    """
    out, i, n = [], 0, len(ts)
    while i < n:
        if ts[i] == "<":
            j = ts.find(">", i)
            if j == -1:
                out.append(ts[i])
                i += 1
            else:
                out.append(ts[i:j + 1])
                i = j + 1
        else:
            out.append(ts[i])
            i += 1
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", required=True, type=Path, nargs="+")
    ap.add_argument("--bags", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--blocks", type=Path,
                    help="blocks.tsv, to mark which called regions coincide "
                         "with a recurrence block")
    ap.add_argument("--members", type=Path)
    ap.add_argument("--model", default="gbrixi/minerva-mlm-8k")
    ap.add_argument("--max-cases", type=int, default=20,
                    help="cases per bag; structure is a per-locus property and "
                         "a handful of loci establishes it")
    ap.add_argument("--threshold", type=float, default=0.6,
                    help="base-pair confidence, Minerva's own default")
    ap.add_argument("--min-region", type=int, default=40,
                    help="noncoding regions shorter than this are not folded; "
                         "matches step 4's --min-len so the two agree on what "
                         "counts as a segment")
    ap.add_argument("--max-tokens", type=int, default=8000)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    import torch
    from minerva.modeling_minerva import MinervaForMaskedLM
    from minerva.rna_structure import call_structures
    import minerva.gene_calling as gc
    import minerva.sequence_utils as su

    def get_tokenizer(mid):
        from transformers import AutoTokenizer
        try:
            return AutoTokenizer.from_pretrained(mid)
        except Exception:
            return AutoTokenizer.from_pretrained(mid, trust_remote_code=True)

    bags = [b for i, b in enumerate(sorted(args.bags.read_text().split()))
            if i % args.nshards == args.shard]
    src = {}
    for root in args.windows:
        if root.is_dir():
            for d in sorted(root.iterdir()):
                if d.is_dir() and (d / f"{d.name}.jsonl").exists():
                    src.setdefault(d.name, d)

    blocks = {}
    if args.blocks and args.blocks.exists():
        for r in csv.DictReader(open(args.blocks), delimiter="\t"):
            if r["decoy"] == "False":
                blocks.setdefault(r["clade"], []).append(r)

    log(f"loading {args.model}")
    tokenizer = get_tokenizer(args.model)
    model = MinervaForMaskedLM.from_pretrained(args.model)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(dev).eval()
    log(f"model on {dev}")

    rows, per_bag = [], []
    for bag in bags:
        d = src.get(bag)
        if d is None:
            log(f"  {bag}: no windows")
            continue
        cases = []
        with open(d / f"{bag}.jsonl") as fh:
            for line in fh:
                cases.append(json.loads(line))
                if len(cases) >= args.max_cases:
                    break
        nreg = npaired = 0
        best = 0.0
        for w in cases:
            seq, ao = w["seq"], w["anchor_offset"]
            inp = gc.build_minerva_input(seq, meta=True,
                                         max_tokens=args.max_tokens)
            enc = su.tokenize_mixed_window(inp["token_string"], tokenizer)
            ids = enc.get("input_ids")
            tokens = split_mixed_tokens(inp["token_string"])
            if ids is None:
                log(f"  unexpected tokenize_mixed_window keys: {list(enc)}")
                return 2
            if enc.get("n_tokens") not in (None, len(tokens)):
                log(f"  token count mismatch: tokenizer {enc['n_tokens']} vs "
                    f"split {len(tokens)}")
                return 2
            # token -> (genome_start, genome_end) for this window. Two traps
            # here, both silent:
            #   * a region's offset from call_structures is in TOKEN space, and
            #     in a mixed locus one amino-acid token spans three bases, so it
            #     is not a window coordinate and must be mapped.
            #   * token_to_genome is indexed by CHARACTER of the token string,
            #     not by token. `<+>`/`<->` are one token but three characters,
            #     so every token after the first strand marker sits at a
            #     different index in the two sequences. Indexing it by token
            #     number silently mismaps 1753 of 2300 nucleotide tokens.
            t2g = inp["token_to_genome"]
            coff, c = [], 0
            for t in tokens:
                coff.append(c)
                c += len(t)
            if c != len(inp["token_string"]) or len(t2g) != c:
                log(f"  token/char mapping inconsistent: chars {c}, "
                    f"token_string {len(inp['token_string'])}, t2g {len(t2g)}")
                return 2
            ids = torch.as_tensor(ids)
            if ids.dim() == 1:
                ids = ids.unsqueeze(0)
            with torch.no_grad():
                out = model(ids.to(dev), output_interactions=True,
                            return_dict=True)
            bp = out.interactions["base_pairing"][0].float().cpu()
            structs = call_structures(bp, tokens, threshold=args.threshold,
                                      min_length=args.min_region)
            for k, st in enumerate(structs):
                pairs = st.pairs or []
                L = len(st.sequence)
                if L < args.min_region:
                    continue
                conf = [p[2] for p in pairs]
                frac = round(2 * len(pairs) / L, 4) if L else 0.0
                t0, t1 = st.offset, st.offset + L - 1
                if t1 >= len(coff):
                    continue
                g0, g1 = t2g[coff[t0]], t2g[coff[t1]]
                if g0[0] < 0 or g1[1] < 0:
                    continue        # region starts or ends on a strand marker
                s0, s1 = g0[0], g1[1]
                hit = ""
                for b in blocks.get(bag, []):
                    bs = ao + int(float(b["median_start"]))
                    be = bs + int(float(b["median_len"]))
                    if min(s1, be) - max(s0, bs) > 0.5 * min(L, be - bs):
                        hit = b["block_id"]
                        break
                nreg += 1
                npaired += 1 if frac >= 0.2 else 0
                best = max(best, frac)
                rows.append({
                    "bag": bag, "anchor_id": w["anchor_id"],
                    "region_index": k,
                    "window_start": s0, "window_end": s1,
                    "anchor_frame_start": s0 - ao,
                    "anchor_frame_end": s1 - ao,
                    "length": L,
                    "minerva_n_pairs": len(pairs),
                    "minerva_paired_fraction": frac,
                    "minerva_mean_pair_conf": (
                        round(sum(conf) / len(conf), 4) if conf else ""),
                    "minerva_max_pair_conf": (
                        round(max(conf), 4) if conf else ""),
                    "minerva_dot_bracket": st.dot_bracket,
                    "coincides_with_block": hit,
                })
        per_bag.append({
            "bag": bag, "n_cases_scored": len(cases),
            "n_noncoding_regions": nreg,
            "n_regions_paired_ge_0.2": npaired,
            "max_paired_fraction": round(best, 4),
            "n_recurrence_blocks": len(blocks.get(bag, [])),
        })
        log(f"  {bag}: {len(cases)} cases, {nreg} regions, "
            f"{npaired} paired, best {best:.3f}")

    sfx = f"_{args.shard:02d}" if args.nshards > 1 else ""
    if rows:
        cols = list(rows[0])
        with open(args.out / f"minerva_regions{sfx}.tsv", "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    if per_bag:
        cols = list(per_bag[0])
        with open(args.out / f"minerva_by_bag{sfx}.tsv", "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in per_bag:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"{len(bags)} bags, {len(rows)} noncoding regions folded -> "
          f"{args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
