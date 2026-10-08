# Blind protein-anchored ncRNA discovery — corrected run

Generated 2026-10-07 from job **26739540**, the first end-to-end run of the
corrected chain: dedup-sequence clade identity → stratified window cap →
quota-constrained within-family decoys → step 4 → benchmark → export → fold →
triage.

Source directories under `/global/scratch/users/kh36969/protein_ncrna/`:
`clades_final_26739540`, `windows_final_26739540`, `blocks_final_26739540`,
`step4_candidates_final_26739540`, `known_retron_benchmark_final_26739540`.

**This supersedes the earlier report**, which is kept as a development log at
[`docs/DEPRECATED_REPORT_raw_anchor_identity.md`](docs/DEPRECATED_REPORT_raw_anchor_identity.md).
Nothing in it is carried forward as a result. In particular the
`GCA_963568795.1 CAUXVC…` GroupII_RT block from that run is **not** treated as a
reproduction target: it belonged to a clade set the broken identity metric
produced.

The discovery path never reads a known RNA. No guide RNA, msr/msd, tracrRNA or
bridge RNA sequence, length or structure model enters anchor selection, clade
decomposition, window extraction, block calling or scoring. RetronDB and Rfam
appear only in `known_retron_benchmark_*/`, after the predictions are frozen,
and only as annotation. Enforced by AST-based tests, not by intent
(`tests/test_leakage.py`).

R-scape resamples its null on every invocation, so covarying-pair counts move by
about one between runs. Each count quoted is the one in the named file.

---

## Summary

> After correcting exact-duplicate protein inflation in clade divergence
> estimates, a stratified cap, and a quota-preserving same-family permutation
> null, the blind protein-anchored pipeline selected 75 clades and called 128
> recurrent noncoding blocks above a 0.80 score floor. No permutation decoy
> reached this score. Among 121 non-Cas3 blocks folded de novo, three had
> adequate R-scape power and significant covariation with zero false-positive
> pairs, including one generic RT clade with a 694-bp upstream block abutting
> the RT ORF.

**Pipeline status: frozen.** Steps 3 and 4 now carry two independent locks — a
quota-constrained same-family null (decoy max 0.719) and adequate-power R-scape
support in three independent blocks. Changing the cap, the reference picker or
the score floor invalidates both and requires a full-chain re-run, so none of
them moves absent a coordinate bug or evidence that the top hits are artifacts.
Work from here is interpretation of the existing candidates, not a search for
higher scores.

---

## 1. Selected clades

`clades_final_26739540/selected_clade_summary.tsv`

Clade identity is now estimated on **distinct protein sequences**: exact
duplicates are collapsed before sampling, a single distinct sequence scores 1.0
by definition, and the pairwise estimate runs over the deduplicated set. The
previous estimate ran over raw anchors, where these clades are mostly duplicates,
so it measured copy number instead of divergence.

| family | rung | clades | nr100 | complete windows | median identity |
|---|---:|---:|---:|---:|---:|
| RT | 0.35 | 33 | 1,411 | 32,175 | 0.564 |
| Cas3 | 0.35 | 6 | 1,413 | 21,226 | 0.708 |
| IS110 | 0.35 | 20 | 1,490 | 20,529 | 0.822 |
| TnpB_Cas12f | 0.70 | 6 | 529 | 12,335 | 0.866 |
| GroupII_RT | 0.50 | 10 | 369 | 5,716 | 0.575 |
| **total** | | **75** | **5,212** | **91,981** | |

75 clades, up from the 43 the broken metric selected. The leading rejection
reason is now `nr100<8` (5,482 clades), with `identity>0.9` second (805) —
previously `identity>0.9` dominated, because duplicate-heavy clades were being
scored as near-identical.

**Do all 75 have enough complete capped windows?** All 75 produced a window set;
19,252 windows were carved, median 197 per clade. 23 clades fall below 50
windows, and that is where the shortfall concentrates:

| family | clades | min | median | max | <50 windows |
|---|---:|---:|---:|---:|---:|
| RT | 33 | 10 | 223 | 500 | 10 |
| IS110 | 20 | 9 | 92 | 500 | 9 |
| GroupII_RT | 10 | 19 | 258 | 500 | 1 |
| Cas3 | 6 | 17 | 165 | 500 | 2 |
| TnpB_Cas12f | 6 | 27 | 500 | 500 | 1 |

9 of 75 clades yielded no real block. Their median window count is 32 against
223.5 for clades that did yield one, and 5 of the 9 have fewer than 35 windows —
so the dominant failure mode is a thin clade, not a failed caller. The two
exceptions (`TnpB_Cas12f GCA_002300265.2` and `RT GCA_053048085.1`, both at the
500-window cap) are genuine negatives.

### Post-cap diagnostics

`windows_final_26739540/cap_diagnostics.tsv`

The cap is only safe if it removes copies rather than diversity. Measured:

| family | clades | capped | median post-cap nr100 | min | few_distinct (<20) | near_identical |
|---|---:|---:|---:|---:|---:|---:|
| RT | 33 | 12 | 25 | 7 | 13 | 0 |
| IS110 | 20 | 5 | 14 | 5 | 13 | 1 |
| GroupII_RT | 10 | 4 | 23 | 9 | 4 | 0 |
| Cas3 | 6 | 1 | 32 | 7 | 2 | 0 |
| TnpB_Cas12f | 6 | 4 | 58 | 6 | 2 | 0 |

26 clades hit the cap; **median non-redundant retention among them is 1.00** —
in 25 of 26 the cap discarded only exact duplicates. The single real loss is
`Cas3 GCA_011214605.1` (792 → 432 distinct proteins, from 20,427 anchors).

34 clades carry fewer than 20 distinct anchor proteins. That is a property of
the corpus, not of the cap, and it is the main reason covariation power is scarce
below. Recording it here is the point: without this table, a clade that arrives
at R-scape with 9 distinct proteins reports "underpowered", which reads like
biology.

## 2. Decoy calibration

`permutation_null.py`, `step4_candidates_final_26739540/step4_family_null.tsv`

225 decoy clades (3 replicates × 75 real clades), matched on family, size, window
length, anchor geometry, ORF density and GC, with only the relatedness between
members removed. Decoys now draw under a per-source quota (`--min-source-clades 4`,
`--max-source-fraction 0.5`), after the uniform draw was found to leave some
decoys up to 71% composed of windows from the single clade they mirrored.

| family | pool | source clades | decoys | degenerate | null strength |
|---|---:|---:|---:|---:|---|
| RT | 8,877 windows | 33 | 99 | 0 | **strong** |
| IS110 | 4,042 | 20 | 60 | 0 | **strong** |
| GroupII_RT | 2,927 | 10 | 30 | 0 | **strong** |
| Cas3 | 1,299 | 6 | 18 | 0 | **strong** |
| TnpB_Cas12f | 2,107 | 6 | 18 | 0 | **strong** |

**Every family has a strong null and zero degenerate decoys.** This is the first
run where that is true — Cas3 previously had too few clades to draw a real null
from, which is what made its FDR 0.556 meaningless.

### Score floors, re-derived

The 0.80 / 0.90 floors in `configs/blind/step4.json` were **re-derived from this
run's decoy distribution**, not carried over. They land on the same two values,
but the justification is new and the old one is withdrawn.

| | n | p50 | p90 | p95 | p99 | max |
|---|---:|---:|---:|---:|---:|---:|
| real blocks | 263 | 0.7900 | 0.9770 | 0.9864 | 0.9920 | 0.9940 |
| decoy blocks | 45 | 0.4387 | 0.6553 | 0.6658 | 0.7194 | **0.7194** |

0.80 sits 0.08 above the highest-scoring decoy block in the entire run, so it is
above the whole observed null rather than at a chosen quantile of it. One caveat
is recorded deliberately: only Cas3 (28 decoy blocks, max 0.7194) and GroupII_RT
(17, max 0.5648) produced any decoy blocks at all. RT, IS110 and TnpB_Cas12f
contributed zero, so they constrain the floor only from above.

## 3. Step-4 candidates

`blocks_final_26739540/discovery_summary.json`, `step4_candidates_final_26739540/step4_candidates.tsv`

| | real | decoy |
|---|---:|---:|
| clades | 75 | 225 |
| clades with a block | 66 | 25 |
| blocks | 263 | 45 |
| median composite | 0.790 | 0.439 |
| p90 composite | 0.977 | 0.655 |
| max composite | 0.994 | 0.719 |

| family | blocks ≥0.80 | decoy blocks ≥0.80 | blocks ≥0.90 | decoy blocks ≥0.90 | null |
|---|---:|---:|---:|---:|---|
| RT | 62 | 0 | 30 | 0 | strong |
| IS110 | 29 | 0 | 18 | 0 | strong |
| GroupII_RT | 27 | 0 | 15 | 0 | strong |
| Cas3 | 7 | 0 | 4 | 0 | strong |
| TnpB_Cas12f | 3 | 0 | 1 | 0 | strong |
| **total** | **128** | **0** | **68** | **0** | |

### How to state this, and how not to

The number this run supports is **not** "FDR = 0". What was measured is:

> Across 225 same-family quota-constrained permutation decoy clades, **no decoy
> block reached 0.80**. At the resolution this empirical null provides, the
> region above 0.80 is decoy-free.

That is an empirical upper bound on the false-discovery rate, not an estimate of
it. 225 decoy clades yielding 45 blocks can show the absence of decoy support at
a threshold; they cannot resolve a true rate below roughly 1/45. A reader who
takes "FDR 0.0" at face value is reading a sample size as a measurement. The
`step4_family_null.tsv` column is named `empirical_fdr_ge_080` and prints 0.0 —
that column is a ratio against an observed count of zero, and should be quoted as
"no observed decoy support", never as a rate.

Within that limit the separation is unambiguous: the two distributions do not
overlap anywhere above 0.72, in any family, including the family whose decoys
score highest.

This answers the question the correction was run to answer: *on the corrected
clade set, does blind recurrence still produce noncoding blocks above the
within-family null?* **Yes, and more cleanly than before** — 128 candidates
against 43-clade-era numbers, with a null that is now trustworthy in all five
families rather than degenerate in one.

### Known-retron benchmark (annotation only, after freezing)

10 of 33 selected RT clades contain a known retron. 9 of those 10 produced at
least one block; 6 produced a block overlapping the annotated ncRNA. Median rank
of the first overlapping block is **1** (4 of 6 at rank 1). Retron-positive RT
blocks score *lower* than unknown RT blocks (median 0.746 vs 0.814), which is the
expected direction: the benchmark is an instrument reading the blind caller, not
a target it was tuned toward.

## 4. Fold and covariation

`step4_candidates_final_26739540/fold_all080/fold_summary.tsv`

121 of the 128 candidates were folded (RNAalifold on soft-ORF-padded blocks,
R-scape two-set test on the core alignment). Verdicts under the frozen rule —
`covariation_supported` requires expected-to-covary ≥3 **and** TP ≥3 **and**
PPV ≥0.80:

| verdict | n |
|---|---:|
| underpowered | 87 |
| covariation_borderline | 29 |
| **covariation_supported** | **3** |
| covariation_test_failed | 2 |

The honest reading of this table is a statement about **power**, not about
structure. Only 3 of 121 blocks reached adequate R-scape power (expected
covarying pairs ≥3); 109 had none at all. **Not one block in the set is
`tested_no_covariation`** — nothing here was tested with adequate power and found
to lack covariation. The 29 borderline blocks all sit below adequate power, so
their 1–5 covarying pairs are suggestive and nothing more.

### The three covariation-supported blocks

| family | clade | side | score | n seq | aln len | mean PID | bp | expected | TP | FP | sens | PPV |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| RT | `GCA_000951815.1 JYCO01000002` | upstream | 0.966 | 200 | 1029 | 0.473 | 34 | 6.3 | **14** | 0 | 41.2 | **100%** |
| GroupII_RT | `GCA_006954725.1 AAGSHD010000039` | upstream | 0.885 | 200 | 1262 | 0.753 | 109 | 15.3 | **23** | 0 | 21.1 | **100%** |
| GroupII_RT | `GCA_002035215.1 MXNX01000032` | upstream | 0.808 | 159 | 853 | 0.949 | 116 | 7.2 | **21** | 0 | 18.1 | **100%** |

All three have **zero false positives** and PPV 100%, and all three retained
100% of their clade's distinct proteins through the window cap (155/155, 104/104,
32/32) — so their power was not manufactured by sampling, and the blocks that
lack power were not capped into it.

Two of the three are GroupII_RT and one is RT; two different families, three
different clades, found blind.

## 5. Locus geometry

`step4_candidates_final_26739540/triage/triage_loci.tsv`, `triage_blocks.tsv`

121 blocks resolve into 110 loci (7 loci carry more than one block). Architecture:
54 upstream, 48 downstream, 6 flanking both sides, 2 split. 20 `nearby_block_group`s
record blocks that abut without being merged into one locus.

The three supported blocks all have clean anchor geometry:

| clade | side | rel. start (q25..q75) | rel. end (q25..q75) | member frac | padded coding frac |
|---|---|---|---|---:|---:|
| RT `GCA_000951815.1` | upstream | −694 (−694..−682) | 0 (−17..0) | 0.966 | 0.178 |
| GroupII_RT `GCA_006954725.1` | upstream | −1089 (−1092..−1016) | −73 (−76..−72) | 0.938 | 0.064 |
| GroupII_RT `GCA_002035215.1` | upstream | −556 (−556..−554) | 0 (0..0) | 0.811 | 0.210 |

Two of them end at **exactly** x = 0 — the anchor's first base — with an
interquartile range of 0 bp across 159 and 483 members. The highest-scoring
blocks in the whole run show the same signature: the two 0.994 GroupII_RT blocks
end at 0 (q25..q75 = 0..0) across 497 members. A block boundary that coincides
with a start codon to the base across hundreds of genomes is the geometry of an
element whose RNA abuts its own ORF, which is what a group II intron with its ORF
masked looks like — and it recurs here independently across GroupII_RT clades.

## 6. What this run does and does not establish

**Does:**
- Blind, protein-anchored recurrence produces 128 noncoding blocks that are
  cleanly separated from a within-family permutation null in all five families,
  with the null now strong everywhere and the floors re-derived from it.
- Three of those blocks have adequate covariation power and show covariation at
  100% PPV with zero false positives.
- The recurrent anchor-abutting geometry survives the correction; it was not an
  artifact of the broken clade set.
- The known-retron benchmark recovers at median rank 1 without ever having been
  consulted during discovery.

**Does not:**
- Establish that the other 118 blocks lack structure. 87 are underpowered and 29
  are borderline; none was tested with adequate power and found negative. The
  limiting resource is distinct sequences per clade (34 of 75 clades have <20),
  not folding or alignment.
- Validate anything biochemically. Every claim here is computational.
- Clear Cas3. It shows no decoy support above the floor on the corrected null,
  but it supplies the highest-scoring decoy blocks in the run, and PF18019
  remains a promiscuous HD domain; `benchmark_only` is retained on those
  grounds, with the old justification withdrawn.
- Measure a false-discovery *rate*. See "How to state this" in section 3: 225
  decoy clades bound the rate from above, they do not estimate it.

---

# Interpretation reports

The pipeline is frozen. What follows reads the frozen candidates; it does not
re-score them. Produced by `scripts/discovery_shortlist.py`,
`scripts/locus_report.py` and `scripts/interpret_locus.py` (the last of which is
UNBLIND and marked so).

## 7. Top RT supported locus — `RT GCA_000951815.1 … B002`

The discovery branch's best result, and the one worth human attention first. The
two GroupII_RT supported blocks are positive-control-like: a structured RNA arm
beside a group II intron RT is the expected answer. This one is a generic RT
clade with no known RNA on it.

`step4_candidates_final_26739540/loci/RT_000951815_B002.*`,
`interpretation/RT_000951815_B002.interpretation.md`

**Clade.** 9,643 anchors capped to 500; **155 distinct proteins before the cap
and 155 after** — the cap cost nothing. Post-cap median identity 0.459, so the
clade is genuinely diverged rather than a duplicate stack. That divergence is
what bought the covariation power.

**Geometry.** 483 members, 96.6% of the clade, **all 483 on the same strand**.
The block spans x = −694 to x = 0: it ends exactly on the RT start codon, with
start IQR 12 bp and end IQR 17 bp across 483 genomes.

The locus table makes a stronger statement than the triage row does. Median
`block_frac_of_run` is **1.000** and median `run_len` is **694** = the block
length: the block is not *inside* an intergenic region, it **is** the entire
noncoding run. Median distance to the upstream ORF is 0 and to the downstream
ORF (the RT itself) is 0. The element occupies the complete gap between the
preceding gene and the RT, exactly, in the median genome.

**Known-RNA overlap.** None. 0 of 483 members overlap any RetronDB ncRNA,
maximum overlap 0 bp. It ranks first in its clade. This clade is not one of the
10 retron-positive RT clades.

**Neighbourhood.** Only 8 of 483 anchors have a scanned marker nearby (HNH,
median 5.1 kb). No cascade genes. On the scanned marker set this does not look
like a CRISPR locus, and there is not enough HNH/TOPRIM/TIR to call it a defence
RT either. Worth stating plainly: only those groups were scanned, so this is
weak evidence of absence, not evidence of a new class.

**Covariation.** 14 covarying pairs, 0 false positives, PPV 100%, and they fall
on **7 distinct helices** rather than piling into one hairpin:

| helix | i | j | pairs | loop span | max power |
|---:|---:|---:|---:|---:|---:|
| 1 | 173 | 664 | 1 | 491 | 0.20 |
| 2 | 179–184 | 213–218 | 4 | 34 | 0.42 |
| 3 | 250 | 261 | 1 | 11 | 0.28 |
| 4 | 298–300 | 358–360 | 2 | 60 | 0.39 |
| 5 | 307–308 | 329–330 | 2 | 22 | 0.20 |
| 6 | 343 | 350 | 1 | 7 | 0.11 |
| 7 | 463–466 | 471–474 | 3 | 8 | 0.28 |

Seven separately supported helices including one long-range pair (173×664) is a
multi-domain fold, not a single stem-loop that happened to covary.

**Reading.** A clade-specific structured noncoding element filling the entire
intergenic gap immediately upstream of a reverse transcriptase, in 96.6% of a
155-distinct-protein clade, on one strand, with covariation support on seven
helices and no overlap with any described retron ncRNA. This is the shape the
project was built to find, and it was found without any RNA prior.

## 8. High-score anchor-abutting blocks

`triage/anchor_abutting.tsv` — score ≥ 0.90, non-Cas3, near edge within 10 bp of
the anchor boundary with edge IQR ≤ 25 bp. **21 blocks**: 10 RT, 9 GroupII_RT,
2 IS110. Typical gap is 0 bp with IQR 0.

A note on the filter, because the obvious column is the wrong one.
`distance_to_anchor` is the offset of the block's *start*. For a downstream
block, 0 means it abuts the stop codon; for an upstream block it is the *far*
edge, and 0 is impossible by construction. Filtering `distance_to_anchor == 0`
selects 7 blocks, **all downstream**, and drops every upstream abutting block —
including both 0.994 GroupII_RT blocks and the supported RT block. So
`anchor_gap_bp` measures whichever edge faces the anchor, with
`anchor_edge_iqr_bp` for its tightness.

The signature is not GroupII-specific: **10 of 21 are plain RT and 2 are IS110**,
across both orientations (12 ending on a start codon, 9 beginning on a stop
codon). A boundary shared with an ORF to the base across hundreds of genomes is
a cis-module geometry independent of covariation power, and most of these 21 are
underpowered rather than negative.

## 9. Discovery triage top 10

`triage/discovery_top10.tsv` — GroupII_RT, Cas3 and known-retron overlaps
removed, leaving a pool of 91 blocks (59 RT, 29 IS110, 3 TnpB_Cas12f); 3 blocks
excluded as known-retron overlaps.

| # | family | score | verdict | power | TP/FP | member frac | side | gap bp | len | coding |
|---:|---|---:|---|---|---|---:|---|---:|---:|---:|
| 1 | RT | 0.966 | **supported** | adequate | 14/0 | 0.966 | up | **0** | 694 | 0.18 |
| 2 | RT | 0.986 | borderline | none | 1/0 | 0.992 | up | 1280 | 178 | 0.46 |
| 3 | IS110 | 0.982 | borderline | low | 3/0 | 0.984 | up | **0** | 70 | 0.68 |
| 4 | RT | 0.971 | borderline | none | 1/0 | 0.992 | up | **0** | 671 | 0.18 |
| 5 | RT | 0.902 | borderline | none | 4/0 | 0.902 | up | 1506 | 319 | 0.30 |
| 6 | IS110 | 0.894 | borderline | none | 1/0 | 0.936 | down | 27 | 237 | 0.13 |
| 7 | RT | 0.886 | borderline | none | 1/0 | 0.886 | up | 348 | 139 | 0.52 |
| 8 | RT | 0.878 | borderline | none | 1/0 | 0.898 | up | 1268 | 386 | 0.28 |
| 9 | RT | 0.856 | borderline | none | 2/2 | 0.856 | down | 178 | 146 | 0.00 |
| 10 | RT | 0.856 | borderline | none | 1/0 | 0.856 | down | 926 | 255 | 0.37 |

Ranked by verdict, then PPV band, then score, member fraction, anchor gap and
coding fraction. The PPV band is not cosmetic: `covariation_borderline` only
asks whether TP > 0, so a block with **TP=1 and FP=118 (PPV 0.84%)** carries the
same label as one with TP=4 and FP=0, and ranking on score alone put the first
above the second. Two such blocks (FP=118 and FP=5) ranked 6th and 4th before
the band was added and are now outside the top 10.

Ranks 2–10 are all underpowered or low-power. They are a queue for deeper
annotation, not nine more results — only #1 has been tested with adequate power.
#4 is the one to look at next: 671 bp, abutting the anchor at 0 bp, 99.2% member
fraction, 18% coding, TP=1/FP=0. Same shape as #1, without the sequence
divergence to prove it.
