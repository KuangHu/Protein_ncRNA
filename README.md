# Protein_ncRNA

> **STATUS: steps 1–6 built and run end to end. The published result numbers
> are DEPRECATED pending a corrected re-run — see the banner in
> [`REPORT.md`](REPORT.md).**
>
> Review found that step 3 estimated clade identity over raw anchors, which in
> these duplicate-heavy families measures copy number rather than divergence.
> Since that statistic gates clade selection, it decided which clades the whole
> pipeline ran on. Corrected (deduplicate on sequence, then sample), the
> selected clade set goes from 43 to 75 and every downstream number changes.
> A corrected end-to-end run is in progress; its output will replace
> `REPORT.md`.
>
> The method, the blind-control discipline and the test suite stand. What is
> withdrawn is the specific set of candidate blocks and the one GroupII_RT
> covariation hit, which must be rediscovered on the corrected clade set before
> it is claimed again.

Blind discovery of protein-associated noncoding RNAs. A protein or locus is the
**only** anchor: the discovery stage uses no known guide RNA, msr/msd, tracrRNA
or bridge RNA sequence, length or structure model. Candidates are then put
through de novo comparative RNA tests and, finally, biochemistry.

Sister projects in `../`: `Fixed_target_editor/` (protein-guided fixed-target
insertion — the same blind-control methodology, applied to DNA targets),
`DL_RNA_guide_edotor_classifer/`, `DL_RNA_guide_edotor_positive_generator/`.

```bash
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python   # py3.11

for t in tests/*.py; do $PY "$t"; done

# 1. corpus accounting -> real_genomes.tsv, the shared entry point
$PY scripts/genomes_db_report.py --out out            # add --verify-gzip on SLURM

# 2. census: anchors + flank feasibility, no window sequence
sbatch sbatch/anchor_census.sbatch                    # 16-shard array
CENSUS_DIR=... sbatch sbatch/census_merge.sbatch      # merge + divergence ladder

# 3. clade decomposition -> selected_clades.tsv, the step-4 entry point
$PY scripts/refilter_census.py --census <census> --out <census>/corrected
sbatch sbatch/neighborhood_scan.sbatch                # Cas3 + RT locus sanity
$PY scripts/clade_decompose.py --census <corrected> --out out/clades \
    --identity --neighborhood <nbr>/neighborhood.tsv

# 4. discovery, on selected clades only -- never family-wide
sbatch sbatch/collect_anchors.sbatch                  # --windows all
sbatch sbatch/extract_windows.sbatch                  # normalized windows + proteins

# 5. the null, then the block caller on real + decoy clades together
$PY scripts/permutation_null.py --windows <win> --out <decoys>   # 3 x 43 decoys
sbatch sbatch/discover_blocks.sbatch
sbatch sbatch/candidate_table.sbatch                  # frozen floors -> candidates

# 6. structure test, triage, report
sbatch sbatch/fold_all080.sbatch                      # export + MAFFT/RNAalifold/R-scape
sbatch sbatch/triage_report.sbatch                    # triage_blocks.tsv, triage_loci.tsv
```

Step 5's score floors are **empirical properties of the decoy null**, not
transferable constants. Change `configs/blind/step4.json` and the null has to be
re-run; the file says so in its own `_rationale` block.

**Never run step 4+ on family-wide windows.** A family-wide IS110 or RT
alignment is biological soup; `selected_clades.tsv` exists to prevent it.

**`genomes_db` is 74% empty download stubs.** Never scan `INDEX.tsv` directly;
use `real_genomes.tsv`. This affects every project reading that mirror.

External binaries are resolved from the directory of the running interpreter
(`protein_ncrna/tools.py`), so an absolute `$PY` path works without activating
the env.

## The workflow

| Step | What | Status |
|---|---|---|
| 1 | Anchor proteins: RT / GroupII_RT / IS110 / Cas3 / Cas9 / Cas12 / Cas13 | **built** — `anchors.py`, `orfs.py` |
| 2 | ±5 kb window per anchor, oriented to anchor strand | **built** — `windows.py` |
| — | Corpus census + divergence ladder | **done** — `census_report.py` |
| 3 | Clade by **protein sequence only** | **built** — `clade_decompose.py`, `neighborhood_scan.py` |
| — | Within-family permutation null: 129 decoy clades | **done** — `permutation_null.py` |
| 4 | Conserved noncoding blocks within a clade | **built** — `discover_blocks.py`, `candidate_table.py` |
| 5 | Candidate export: core + soft-ORF-padded FASTAs | **built** — `export_blocks.py` |
| 6 | De novo comparative RNA test: MAFFT → RNAalifold → R-scape | **built** — `fold_blocks.py` |
| — | Locus triage and per-locus reports | **built** — `triage_report.py`, `locus_report.py` |
| 7 | Exclude coding / repeat / UTR artifact | not started |
| 8 | RNA-seq / RACE boundaries | wet lab |
| 9 | Deletion, catalytic-dead, trans rescue, stem compensation | wet lab |

Step 5 is a **pluggable proposer**: local repeat-block and paired-block
detectors first; Minerva (`garykbrixi/minerva`) can be added later as another
scorer. Minerva proposes candidates — it is never the final evidence.

### Why a ladder, not one threshold

`nr90` showed the families are not all duplicates, but RNA discovery needs
clades in a narrower band:

| protein identity | expected behaviour |
|---|---|
| > 95% | RNA nearly identical; R-scape has no substitutions to score |
| 70–90% | usually good for a tight element subfamily |
| 35–70% | may keep the scaffold, may lose alignability |
| < 35% | normally too remote unless the RNA is severely constrained |

The right rung differs per family, so `clade_decompose.py` clusters at
90/70/50/35% and `recommended_thresholds.json` picks one rung per family.
Clades at different rungs are nested, so step 4 must consume exactly one rung
or it will process the same locus repeatedly.

### Anchor priority after the census

| priority | family | why |
|---|---|---|
| 1 | **TnpB** | known ωRNA immediately downstream of *tnpB*; anchor-to-RNA geometry is the clearest available |
| 2 | **IS110** | bridge RNA is the most on-target guide-RNA control, and the most anchors |
| 3 | **GroupII_RT** | highest information density (62× redundancy); RT + structured RNA |
| 4 | RT | must first be split into retron / Abi-like / group-II-like / host |
| 5 | Cas3 | *provisional* — candidate anchor family, not a control, until cascade genes are confirmed |
| drop | Cas9, Cas13 | this corpus cannot support them |

RT subtypes are labelled only coarsely (`unknown_compact`, `groupII`,
`candidate_defence`, `long_maturase_like`) — enough to stop incompatible loci
being pooled. **"Retron" is never claimed from RT homology alone**, since
generic RT covers retrons, Abi-like defence RTs, DGRs, group II remnants and
host enzymes alike.

## Keeping the blind

The three ground-truth systems are blind positive controls, so the discovery
stage must not be told the answer:

| forbidden in discovery | allowed as anchor |
|---|---|
| Rfam retron RNA CM | RT protein HMM |
| known bridge RNA coordinates | IS110 protein anchor |
| tracrRNA annotation | Cas9 protein anchor |
| msr/msd length or fold rule | ±5 kb raw sequence |
| crRNA repeat consensus | de novo repeat block |

Enforced mechanically, not by intention:

- `configs/anchors.json` holds protein models only, and `load_config` raises on
  any accession beginning `RF`.
- `tests/test_blind.py` fails if a discovery module imports `unblind/` or
  hardcodes an Rfam accession.
- Known RNA coordinates live in `unblind/` and `data/known_retrons/`. Nothing
  under `protein_ncrna/` may read either. Unblinding
  (`blind prediction ∩ known RNA coordinate`) is a separate script, run only
  after predictions are frozen.
- `tests/test_leakage.py` parses every script's AST and fails if discovery code
  references RNA structure tooling. Step 6 is exempted by an explicit one-entry
  allowlist (`STEP6_SCRIPTS`) with a stale-entry guard — an allowlist, not a
  convenience exemption.
- The known-retron benchmark is an **instrument, not a target**. No parameter is
  ever tuned until a known positive comes back; doing so would convert blind
  discovery into a retron detector and the failure would be silent, because the
  benchmark would then always look good.

## Why step 2 orients to the anchor strand

Contig orientation is arbitrary with respect to biology. Without normalization,
"200 bp upstream of the anchor" means a different thing in each genome and a
genuinely conserved block is smeared across both strands — the sister project
measured roughly halved information content from exactly this. After
`extract_window`, left-of-anchor means the same thing in every record, so a
block recurring at a consistent offset is evidence rather than an artifact.
`tests/test_windows.py` asserts a locus encoded on either strand yields an
identical oriented window.

## Candidate score

Evidence lines, weakest to strongest. No single one is sufficient:

| evidence | what it buys |
|---|---|
| co-occurs within a protein clade | not a random intergenic hairpin |
| same relative locus position | inherited with the element |
| paired / repeat block | the DNA has an RNA-like dependency |
| **covariation** | the secondary structure is supported by evolution |
| TSS + 3′ end | there really is an independent RNA |
| no conserved smORF | not a tiny protein |
| deletion + trans rescue | the function is in the RNA molecule |

## Ground truth

Controls are in **two tiers**, because the original trio leaned too hard on this
corpus — Cas9 returned zero hits in the pilot:

| tier | purpose | depends on genomes_db |
|---|---|---|
| **in-corpus** | blind discovery; chosen *after* the census, needs ≥ 20 diverse loci | yes |
| **external curated** | prove the algorithm recovers known ncRNA, on separately downloaded benchmark loci | no |

| test | anchor | should recover | standing |
|---|---|---|---|
| IS110 | IS110/IS1111 transposase | bridge RNA | retained — diversity lives in the element |
| retron | RT ORF | msr/msd region | retained, not promised — awaiting census |
| group II intron | maturase (GIIM) | the intron RNA itself | **added** — natural RT + RNA test |
| Cas3 / Cas12 / Cas13 | type I / V / VI effectors | CRISPR array ± scaffold | **added** — both Cas12 and Cas3 hit where Cas9 did not |
| Cas9 | Cas9 protein | array + tracrRNA | still scanned; **dropped as an in-corpus control** |

A curated control failing means the method is wrong. An in-corpus control
failing may only mean the corpus is thin — which is why they are kept apart.

## Readouts for a novel locus

Minimum set: WT locus, RNA deletion, protein catalytic-dead, RNA trans rescue,
stem-disrupt, stem-compensatory.

| suspect | catalytic-dead | first readout |
|---|---|---|
| IS110-like | mutate Tnp catalytic residue | donor × target recombination |
| retron-like | mutate RT YADD | msDNA / cDNA hairpin product |
| Cas9-like | dCas9 / nuclease-dead | target cleavage or repression |
| unknown RT | RT-dead | small nucleic-acid product + sequencing |

Deletion inactivates → the candidate is necessary. Trans rescue restores → the
function is an RNA molecule, not just a DNA element. Stem-disrupt inactivates
and the compensatory mutation restores → the function depends on the predicted
structure.

## Outputs

`collect_anchors.py` writes `windows.jsonl` (one oriented window per anchor,
with `anchor_offset`, flank sizes and truncation flags), `anchors.tsv` (flat
table for inspection) and `summary.json` (per-family counts, failures).
