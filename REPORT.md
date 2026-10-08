> # ⚠ DEPRECATED — raw-anchor identity metric
>
> **Every number in this file is invalidated and is retained as a development
> log, not as a result.**
>
> Step 3 estimated each clade's median pairwise identity over *raw anchors*.
> These clades are mostly duplicates — one GroupII_RT clade has 799 members and
> 29 distinct proteins — so that statistic measured copy number, not divergence,
> and converged on 1.0 for duplicate-heavy clades regardless of how diverged
> they were. The same clade measured **0.451** from an ordered sample and
> **1.000** from a random one. Median identity gates clade selection, so a
> broken estimate decided which clades the entire pipeline ran on.
>
> | superseded output | status |
> |---|---|
> | 43 selected clades | invalidated — corrected run selects 75 |
> | 24 non-Cas3 candidates | invalidated — derived from the wrong clade set |
> | the GroupII_RT block with 19 covarying pairs | **not trusted until rediscovered** |
> | the decoy FDR table | invalidated — real clade set changed |
> | Cas3 `benchmark_only` at FDR 0.556 | invalidated — its null was degenerate or absent |
>
> Also fixed since: a minus-strand member-trim sign error, a uniform decoy draw
> that left decoys up to 71% drawn from a single real clade, an encounter-order
> window cap, and a step-4 worker failure path that exited 0.
>
> The corrected chain is dedup-sequence identity + quota decoys + stratified
> cap. Its results will replace this file. Do not cite anything below.

---

# Blind protein-anchored ncRNA discovery — results to date

Generated 2026-10-07 from the post-strand-fix run (`blocks_strandfix_26738251`,
`step4_candidates_strandfix_26738251`). All numbers below come from files under
`/global/scratch/users/kh36969/protein_ncrna/`; each section names its source.

R-scape resamples its null on every invocation, so covarying-pair counts move by
about one between runs. Where a count is quoted, it is the one in the named file,
not a constant.

The discovery path never reads a known RNA. No guide RNA, msr/msd, tracrRNA or
bridge RNA sequence, length or structure model enters anchor selection, clade
decomposition, window extraction, block calling or scoring. RetronDB and Rfam
appear only in `known_retron_benchmark_*/`, after the predictions are frozen,
and only as annotation. This is enforced by AST-based tests, not by intent
(`tests/test_leakage.py`, 21 tests passing).

---

## 1. Anchors and selected clades

`census_26708202/`, `clades_v1/selected_clades.tsv`, `windows_26714415/`

343,125 protein anchors were collected by Pfam HMM search across the genome
corpus and decomposed into clades at a per-family identity threshold chosen by
`recommend_thresholds()` — one rung per family, picked for the largest number of
selected clades, with ties broken toward the *less* identical rung so a tight cut
cannot silently discard diversity.

| family | selected clades |
|---|---:|
| RT | 22 |
| IS110 | 12 |
| TnpB_Cas12f | 4 |
| GroupII_RT | 3 |
| Cas3 | 2 |
| **total** | **43** |

Each selected clade yields one normalized window set in the anchor frame:
`x = window_position − anchor_offset`, so x=0 is the anchor's first base and
negative coordinates are upstream. Windows average ~86% coding.

## 2. Decoy calibration

`permutation_null.py`, `blocks_strandfix_26738251/`, `step4_candidates_strandfix_26738251/step4_family_null.tsv`

The null is a **within-family permutation**: 129 decoy clades (3 replicates ×
43 real clades), matched on family, size, window length, anchor geometry, ORF
density and GC, with only the relatedness between members removed. Decoys go
through the identical block caller. Every operating point in
`configs/blind/step4.json` was set from this null and from nothing else.

| family | real ≥0.80 | decoy (scaled) ≥0.80 | FDR | real ≥0.90 | decoy ≥0.90 | benchmark_only |
|---|---:|---:|---:|---:|---:|:--|
| RT | 15 | 0.0 | 0.000 | 10 | 0.0 | |
| IS110 | 4 | 0.0 | 0.000 | 1 | 0.0 | |
| GroupII_RT | 4 | 0.0 | 0.000 | 3 | 0.0 | |
| TnpB_Cas12f | 1 | 0.0 | 0.000 | 1 | 0.0 | |
| **Cas3** | 3 | 1.67 | **0.556** | 2 | 0.0 | **yes** |

Cas3 is the one family whose decoys clear the score floor, and PF18019 is a
promiscuous HD domain, so Cas3 is flagged `benchmark_only` and excluded from the
candidate set. Across the other four families **no decoy clade produced a block
at or above 0.80.**

Frozen defaults: `diverse_kmer` references, k=3, `min_recurrence` 0.5,
`score_floor` 0.80, `high_confidence_floor` 0.90.

## 3. The 24 non-Cas3 candidates

`step4_candidates_strandfix_26738251/step4_candidates.tsv`, `triage/triage_blocks.tsv`, `triage/triage_loci.tsv`

24 blocks ≥0.80 across 13 clades, grouped into **21 loci**. Grouping is
anchor-geometric: blocks of one clade merge when both lie within 500 bp of the
anchor ORF or abut opposite sides of it. Block-to-block adjacency is reported
separately as `nearby_block_group_id` and deliberately does not merge loci — two
blocks 1 kb from the anchor may be one distal RNA or may be another gene's UTR
and terminator, and calling them a protein-associated locus would assert the
first without evidence.

| architecture | loci |
|---|---:|
| upstream of the anchor ORF | 10 |
| downstream | 7 |
| **flanking both sides** | **3** |
| split (members disagree on side) | 1 |

Only one locus is genuinely `split` (`RT_C10_L1`). Three earlier split calls were
an artifact of testing each member against the clade's *median* anchor length
instead of its own; fixed, and the rule is documented in `side_of()`.

## 4. GroupII_RT positive-control recovery

`triage_loci.tsv`, `locus_groupII/`, `locus_groupII_c01/`

**Both** GroupII_RT clades that produced blocks produced the same architecture,
independently: a long structured noncoding block ending exactly at the RT start
codon and a short one beginning exactly at the stop codon, both on the anchor's
strand.

| | C02 (GCA_963568795.1) | C01 (GCA_003903225.1) |
|---|---|---|
| anchor ORF | 1,389 bp | 1,818 bp |
| upstream arm | −449 → 0, 449 bp, 493 members | −119 → 0, 119 bp, 483 members |
| downstream arm | +1389 → +1465, 76 bp, 478 members | +1822 → +1893, 71 bp, 429 members |
| coordinate scatter | IQR ±1 bp | IQR ±1 bp |
| strand agreement | 493/493 and 478/478 `+` | 483/483 and 429/429 `+` |
| covariation | **19 TP / 4.1 expected / PPV 100 / FP 0** | 3 TP, power `none` |

That is the architecture of a group II intron with its ORF masked out: the RT
sits inside domain IV, domains I–III lie 5′ of it and domains V–VI lie 3′. The
pipeline was not told this. It recovered both arms of one structured mobile RNA
from the protein alone, twice, in unrelated clades.

Caveats kept on the record: the 5′ arm is 449 bp where DI–DIII is typically
900–1,500 bp, truncated by a prodigal ORF call upstream rather than by the
recurrence signal; and nothing here *identifies* the element — no intron model
was consulted, and confirming identity is a benchmark-branch question.

For reference, the frozen retron benchmark (`known_retron_benchmark_26714577/`,
annotation only): of 7 retron-positive RT clades, 6 produced a block, 4 had a
block overlapping the known ncRNA, **median rank 1**, three at rank 1.

## 5. Top unknown RT candidates

`triage_blocks.tsv`

| locus | score | side | coordinates | members | identity | verdict |
|---|---:|---|---|---:|---:|---|
| RT_C09_L1 | 0.946 | **flanks both sides** | −114→0 and +1101→+1659 | 473 / 468 | 0.90 / 0.72 | borderline |
| RT_C09_L2 | 0.938 | upstream | −3755 → −2361 | 469 | 0.88 | underpowered |
| RT_C11_L1 | 0.834 | upstream | −1166 → −1033 | 417 | 0.94 | underpowered |
| RT_C08_L1 | 0.812 | upstream | −2403 → −1907 | 218 | 0.57 | borderline |
| RT_C07_L1 | 0.808 | upstream | −365 → 0 | 19 | 0.60 | underpowered |

RT_C09_L1 is the most interesting unknown: a third clade showing the bipartite
ORF-flanking architecture, 473 of ~490 members, coordinates abutting both ends of
the anchor ORF. Its covariation is thin (1 TP against 2 FP) and should not be
cited as support.

RT_C12 contributes 7 of the 24 blocks from a single 40-window clade at 98–100%
member identity. Those are high-recurrence by construction and carry no
structural information; they should be read as one clade's worth of evidence, not
seven.

## 6. Most candidates are covariation-underpowered, and that is not a negative

`step4_candidates_strandfix_26738251/fold_all080/fold_summary.tsv`

24 blocks folded de novo (MAFFT → RNAalifold → R-scape) on their own members
only. No curated family or covariance model is read at this stage.

| verdict | n |
|---|---:|
| covariation_supported | **1** |
| covariation_borderline | 5 |
| underpowered | 18 |
| tested_no_covariation | **0** |

`supported` requires R-scape expected-to-covary ≥3 **and** ≥3 covarying pairs
**and** PPV ≥0.80. `borderline` is any real covarying pair that misses those.
`tested_no_covariation` — the only genuinely negative verdict — is reserved for
alignments that had the power to show covariation and did not, and **nothing in
this set qualifies**.

The reason is power, measured by R-scape itself:

| power (expected pairs to covary) | blocks |
|---|---:|
| adequate (≥3) | 1 |
| low (1–3) | 1 |
| none (<1) | 20 |
| none — alignment has no substitutions at all | 2 |

**One block out of 24 had the sequence diversity for the test to mean much.** A
null result here is a statement about the alignment, not about the RNA. Two
blocks are 100% identical across all members; R-scape omits its power section
entirely for those.

---

## What can and cannot be claimed

Supportable:

> In a blind protein-anchored ncRNA discovery benchmark, the pipeline recovered
> 24 non-Cas3 candidate blocks above an empirical within-family permutation
> floor, forming 21 loci. Two independent GroupII_RT clades yielded the expected
> bipartite architecture — structured noncoding blocks immediately flanking the
> RT ORF on both sides — and one of them carries strong R-scape covariation
> support (19 covarying pairs against 4.1 expected, PPV 100%). Most remaining
> candidates are too sequence-similar for covariation to be informative, so they
> are reported as underpowered rather than negative.

Not supportable: "5 loci have covariation support"; "10 blocks failed"; any claim
that a novel RNA has been identified. The demonstrated result is that the
pipeline can recover both arms of a structured mobile RNA from a protein anchor
after the ORF is masked — a method result, not a discovery.

## Next

Manual inspection of the top RT, IS110 and TnpB loci — sequences, folds and
member coordinate plots (`locus_*/`*`.relative_coordinate_plot.svg`). Further
detector development has diminishing returns at this point.
