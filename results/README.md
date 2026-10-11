# Results from the two external RNA models

Everything Minerva and RNAMunin produced on this corpus, plus the tables the
benchmarks read alongside them and the benchmark output itself. The pipeline's
own working tree lives on scratch and is not in the repo; this directory is
the part the conclusions in `docs/` rest on.

Regenerate the whole directory with `sbatch sbatch/publish_results.sbatch`.
The score tables are copied; the benchmark tables are re-run, so what is
committed is a run's output and not a transcribed number.

## `minerva/` -- RNAfold-free folding model, per non-coding region

`minerva_paired_fraction` is the score; a bag's score is the **mean** over its
regions, never the max (see `docs/minerva_aggregation.md` -- the max costs
0.16 of AUC).

| file | rows | what |
|---|---|---|
| `bench_regions.tsv` | 3,161 | the six-family benchmark: IS110, IS200/605, group II vs IS1, IS3, IS4/5/10 |
| `bench_by_bag.tsv` | 179 | the same, one row per bag |
| `broad_regions.tsv` | 1,167 | the other six RNA-free families: ISL3, IS256, IS66, Mu, Tn3, IS630 |
| `broad_by_bag.tsv` | 55 | |
| `scan_regions.tsv` | 2,362 | the open scan over every one-or-two-protein element, labelled or not |
| `scan_by_bag.tsv` | 125 | |

Shards are concatenated under one header; the `minerva_*_NN.tsv` split is a
SLURM array artefact and carries no meaning.

## `rnamunin/` -- Rfam-trained sequence CNN, per base

| file | what |
|---|---|
| `nc_munin_scores.tsv.gz` | 18,798 records: every candidate and cluster member stage 1 emitted, scored on both strands |
| `selftest_scores.tsv` | the torch port against an independent numpy forward pass (agreement 2e-6) |

A record is one sequence; `mean_*` is the per-base mean, `best50_*` the best
50-nt window, `_fwd` and `_rc` the two strands. The published weights were
rebuilt in torch because TensorFlow will not import on this cluster.

## `derived/` -- what the benchmarks read next to the models

| file | what |
|---|---|
| `ir_flag.tsv` | longest inverted-repeat arm per candidate. A diagnostic, never a filter in the pipeline |
| `orientation.tsv` | each candidate's placements in the anchor-oriented windows: `frac_plus`, the side of the transposase, whether it had to be flipped |
| `scan_1_2_protein.tsv` | the open scan, 353 bags ranked by the blind score |
| `scan_shortlist.tsv` | its top slice |

## `benchmarks/` -- the numbers in `docs/`

| file | script | question |
|---|---|---|
| `separation_bench.txt` | `separation_bench.py` | every axis against the six-family negative |
| `munin_diagnostics.txt` | `munin_diagnostics.py` | is the RNAMunin score length or GC? |
| `broad_is_bench.txt` | `broad_is_bench.py` | length x Minerva against all twelve families |
| `rnamunin_broad_bench.txt` | `rnamunin_broad_bench.py` | RNAMunin against all twelve, 18 aggregation cells |
| `final_separation_bench.txt` | `final_separation_bench.py` | length x Minerva x strand, identical bag sets |

The family labels come from Pfam domains and are used only to **read** the
answer. Nothing in the pipeline that produces these scores sees a label, a
known guide RNA, or a known msr/msd, tracrRNA or bridge RNA sequence.
