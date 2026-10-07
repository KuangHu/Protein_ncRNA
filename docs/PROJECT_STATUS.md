# Protein_ncRNA — project status

**Started 2026-10-06.** Steps 1–2 are built and validated. The current task is
a **corpus census**, not discovery: steps 3–9 stay unwritten until the census
says which anchor families have enough diversity to support a comparative RNA
test. Read the blockers below first.

Decision order, deliberately not negotiable:

```
fix corpus accounting
  → full anchor census over all real genomes
    → pick clades with real diversity
      → install R-scape + Infernal
        → steps 3-9
```

Running a partial discovery job now would only demonstrate the corpus bias more
expensively, and still would not say whether retron or CRISPR controls are
feasible.

## The design

Protein-anchor → blind ncRNA → biochemical validation. The governing rule:

> Only protein / locus is used as an anchor. The discovery stage uses no known
> guide RNA, msr/msd, tracrRNA or bridge RNA sequence, length or structure model.

Nine steps: (1) anchor proteins → (2) ±5 kb oriented windows → (3) protein-only
clades → (4) conserved noncoding blocks → (5) candidate proposal → (6) de novo
comparative RNA test → (7) coding/repeat/UTR exclusion → (8) RNA-seq/RACE
boundaries → (9) deletion, catalytic-dead, trans rescue, stem compensation.

Blind positive controls, run before any novel locus is believed: IS110 → bridge
RNA, retron RT → msr/msd, Cas9 → CRISPR array + tracrRNA. Unblinding is
`blind prediction ∩ known RNA coordinate`, scored on overlap, rank and clade
coverage.

## Blocker 1 — corpus accounting was wrong, and is now fixed

`INDEX.tsv` has 305,483 rows, but **74% are 20-byte gzip stubs** from failed
downloads. The real count is **79,554 genomes** (files ≥ 100 KB).

An earlier figure of 47,899 was itself too low: it counted only genomes inside
the two IS110-biased species views. **31,655 real genomes sit outside any view**
— roughly 40% of the corpus, and plausibly its most taxonomically varied part,
since the views were built by selecting for IS110.

| species (view label) | indexed | real | % real |
|---|---|---|---|
| *(not in any view)* | 67,231 | 31,655 | 47.1% |
| salmonella_enterica | 116,106 | 22,753 | 19.6% |
| klebsiella_pneumoniae | 82,423 | 13,433 | 16.3% |
| enterococcus_faecium | 29,998 | 7,360 | 24.5% |
| enterobacter_hormaechei | 6,910 | 2,143 | 31.0% |
| burkholderia_pseudomallei | 1,579 | 1,579 | 100% |
| mycobacterium_tuberculosis | 1,040 | 444 | 42.7% |
| leptospira_interrogans | 142 | 137 | 96.5% |
| streptococcus_suis | 47 | 47 | 100% |
| yersinia_enterocolitica | 7 | 3 | 42.9% |

`scripts/genomes_db_report.py` writes **`real_genomes.tsv`**, which is the
intended shared entry point for every pipeline reading this mirror. Do not scan
`INDEX.tsv` directly.

### This contaminates other projects

Any pipeline sampling `INDEX.tsv` naively gets a sample that is ~74% dead files:

| effect | symptom |
|---|---|
| most sampled genomes empty | effective n far below nominal n |
| species depth overstated | thousands claimed, tens real |
| downstream hit rate depressed | wrong denominator |
| no clustering or covariation power | same few real files redrawn |

`DL_RNA_guide_edotor_classifer` and the `fna_ins_discovery` pipelines read this
mirror and should be checked against `real_genomes.tsv`.

### Genomes outside the views had no species label

40% of real genomes are in no view, so they would be counted as one giant
`unknown` species. `collect_anchors.py` now parses the organism from each
assembly's own FASTA header and normalizes it to the same slug the views use,
so a view-labelled and a header-labelled genome of one species collapse to a
single entry instead of inflating the diversity count.

### Diversity, not count, is the question

Covariation needs *diverged* sequences. 13,000 IS110 hits across 13,000
near-identical K. pneumoniae strains is one data point in a large hat, and
R-scape would find nothing in it — silently, for a true RNA as readily as a
false one. `scripts/census_report.py` therefore dereplicates each family at
100 / 95 / 90 / 70% identity: a family whose count barely falls across that
ladder is genuinely diverse; one that collapses is a strain artifact. Families
with ≥ 20 nonredundant sequences at 90% identity are flagged
`covariation_feasible`.

## Blind controls, restructured into two tiers

The original IS110 / retron / Cas9 trio leaned too hard on this corpus. Cas9
returned **zero hits** in the pilot — these enterobacteria largely lack it — and
forcing that run would yield an uninformative negative. Cas9 is not a bad
control; this mirror is just not where Cas9 lives.

| tier | purpose | depends on genomes_db |
|---|---|---|
| **in-corpus positive controls** | blind discovery, chosen *after* the census, requiring ≥ 20 diverse loci | yes |
| **external curated controls** | prove the algorithm recovers known ncRNA, on separately downloaded benchmark loci | no |

Keeping them apart is cleaner: `genomes_db` does discovery, the curated
benchmark proves sensitivity. A curated control failing then means the method is
wrong; an in-corpus control failing may only mean the corpus is thin.

Anchor families now scanned, and their standing:

| family | status |
|---|---|
| `IS110` | retained — plausibly feasible; diversity lives in the element, not the host |
| `RT` (retron falls here) | retained, not promised — awaiting census |
| `GroupII_RT` | **added** — maturase whose partner RNA is the intron itself; a natural RT + RNA blind test |
| `Cas3` | **added** — type I, common in enterobacteria, likely the best-populated CRISPR anchor here |
| `Cas12`, `Cas13` | **added** — broader CRISPR sampling than Cas9 alone |
| `Cas9` | still scanned (a count is informative), but **dropped as an in-corpus control** |

Early signal from a 12-genome smoke test: Cas12 and Cas3 both produced hits
where Cas9 produced none, so broadening the anchor set was worth it.

Group II maturases hit `GIIM` **and** `RVT_1`, so family assignment is
priority-first, score-second — otherwise one biological family would scatter
between `GroupII_RT` and `RT` depending on which domain scored higher in each
sequence. Every model and family that hit stays on the record
(`all_models`, `all_families`, `ambiguous`).

## Blocker 2 — two required tools are not installed, and that is fine for now

| tool | step | status |
|---|---|---|
| R-scape | 6, covariation | **missing**, bioconda-installable |
| Infernal (`cmbuild`/`cmsearch`) | 6, CM build and search | **missing**, bioconda-installable |
| Minerva | 5, contact scoring | not installed; treated as an optional pluggable scorer |

Present and used: HMMER3, Pfam-A (`/global/scratch/users/kh36969/pfam_db`),
prodigal, mmseqs, mafft, ViennaRNA (`RNAalifold`, `RNAfold`), blastn, easel.

Without R-scape, step 6 degenerates to RNAalifold alone. RNAalifold will return
a confident consensus structure for almost any AT-rich bacterial intergenic
alignment, so on its own it cannot separate "conserved noncoding block" from
"conserved because the flanking protein is conserved". Covariation is the step
that makes the blind controls interpretable; it is not optional polish.

**They are deliberately not installed yet.** The census needs only prodigal,
hmmsearch, mmseqs and python. Installing the RNA tools first invites the
embarrassing outcome of a working R-scape with no clade that has enough
informative substitutions to feed it. Install them once the census names the
clades worth the effort — see `scripts/install_env.md`.

Minerva (`garykbrixi/minerva`) is deliberately not a dependency. Step 5 is a
pluggable proposer interface; the first implementation will be local
repeat-block and paired-block detectors, and Minerva can be added as another
scorer later. Per the design, Minerva proposes candidates — it is never the
final evidence.

## Keeping the blind

`configs/anchors.json` carries protein models only, and `load_config` refuses an
accession beginning `RF` so an Rfam model cannot be added by accident. Known RNA
coordinates belong in `unblind/`, which **no module under `protein_ncrna/` may
import**. The unblinding comparison is a separate script run after predictions
are frozen.

## Census results — 2026-10-06, job 26708202

**343,124 anchors over all 79,554 real genomes.** 16-shard array, ~30 min per
shard, 0 empty stubs, 23 gunzip failures (0.03%). The independent gzip
verification found 19 truncated files, agreeing with those 23 to within noise;
the verified manifest is at `corpus_verified/real_genomes.tsv`.

| family | anchors | genomes | species | nr90 | redundancy | complete window | covariation? |
|---|---|---|---|---|---|---|---|
| IS110 | 138,897 | 45,816 | 20 | **1,490** | 93× | 0.76 | **yes** |
| RT | 108,348 | 49,788 | 21 | **1,100** | 99× | 0.84 | **yes** |
| TnpB (was "Cas12") | 51,684 | 25,876 | 14 | **240** | 215× | 0.79 | **yes** |
| Cas3 | 31,869 | 31,496 | 10 | **73** | 437× | 0.85 | marginal |
| GroupII_RT | 12,090 | 10,397 | 12 | **194** | **62×** | 0.71 | **yes** |
| Cas9 | 236 | 198 | 2 | 17 | 14× | 0.79 | **no** |
| Cas13 | 0 | — | — | — | — | — | absent |

Five families clear the diversity bar. The corpus is usable after all — the
earlier pessimism was based on a pilot that saw only the IS110-biased views.

### Reading these numbers

**`nr90`, not `anchors`, is the real n.** IS110 has 138,897 hits but 1,490
nonredundant sequences at 90% identity — a 93× redundancy that is exactly the
conspecific-strain inflation predicted. 1,490 is still ample.

**GroupII_RT is the best anchor per unit effort**: 62× redundancy, the lowest of
any family, so its 12,090 anchors carry more independent information than
Cas3's 31,869 (437×). It was added this round and is already the cleanest
RT + RNA blind test available here.

**Cas3 is marginal.** 73 nr90 clears the ≥20 threshold but falls to 57 at 70%
identity, and 437× redundancy means it is essentially one locus copied across
enterobacteria. Treat a Cas3 result as weak evidence.

**Cas9 is settled: not usable here.** Not zero as the pilot suggested — 236
anchors exist — but 17 nr90 across 2 species cannot support covariation.

**These are family-level numbers, not clade-level.** Covariation needs diverged
sequences *within one alignable clade*. 1,490 nonredundant IS110s are spread
across many IS110 subfamilies that may not be mutually alignable. The per-clade
counts are what decide whether R-scape has power, and producing them is step 3.
`covariation_feasible` is a green light to attempt clade decomposition, not a
promise that step 6 will work.

## Two false-anchor classes the census exposed

**"Cas12" was TnpB.** 51,683 of 51,684 hits came from one model,
`Cas12f1-like_TNB` (PF07282): median 402 aa, ~2 copies per genome. Cas12a is
~1,300 aa and single-copy. PF07282 is the **TnpB** TNB domain — the IS200/IS605
transposon protein that Cas12f descends from. Renamed `TnpB_Cas12f`, with true
Cas12 models kept in a separate `Cas12` family (1 census hit).

This is a promotion, not a demotion. TnpB is ωRNA-guided, and the ωRNA is
encoded immediately downstream of *tnpB* — element-encoded guide RNA beside its
protein, which is precisely this project's target. With 240 nr90 across 14
species it is a strong in-corpus control candidate.

**RVT_3 was admitting a host gene.** 11.5% of RT anchors were ≤160 aa, and
10,784 of them are *the same 138 aa sequence*; the top 5 distinct sequences
cover 11,804 of 12,429. A sequence identical across thousands of genomes is a
host housekeeping gene (the sequence is RNase H-like), not a mobile element, and
138 aa is far too short for a functional RT. Families now carry their own
`min_aa_len`, applied after assignment — RT ≥ 250, TnpB ≥ 250, Cas12 ≥ 800 —
since the right floor differs per family. Dereplication had already absorbed
most of this, so the `nr90` figures above stand.

**Still unverified:** `Cas3_HD` (PF18019) supplies 31,386 of 31,869 Cas3
anchors, and the HD clan is promiscuous. Median length 888 aa matches real Cas3
(~900 aa) and ~1 copy/genome fits a single chromosomal CRISPR locus, so it is
probably sound — but confirm by checking co-location with other *cas* genes
before trusting a Cas3 result.

## Steps 1–2: built and validated

`scripts/collect_anchors.py` runs genomes → prodigal ORFs → anchor HMM search →
oriented ±5 kb windows, writing `windows.jsonl`, `anchors.tsv`, `summary.json`.

Validated on a species-stratified pilot drawn from `genomes_db`. Orientation was
checked against real data, not just the synthetic test: slice each window at
`anchor_offset` for `anchor_len` and the result must be a clean ORF. Across 109
pilot windows, **every complete anchor had a canonical start and stop codon and
was in frame**; the 13 exceptions were all prodigal partial genes running off a
contig edge, which is now carried explicitly as `anchor_partial` rather than
rediscovered downstream. 45/109 anchors were on the minus strand, so orientation
was doing real work.

Observed anchor yield is what the corpus predicts: IS110 and RT are abundant,
**Cas9 has not appeared at all** — these enterobacteria largely lack it. The
Cas9 blind control cannot be run on this mirror.

Performance: prodigal dominates. `-p single` with a `-p meta` fallback is ~3×
faster than `meta` alone (19 s vs 53 s on a 7 Mb genome) at equivalent gene
counts. Budget ~19 s/genome/core; 1,800 genomes on 32 cores is ~20 min.

Run it on SLURM — `sbatch/collect_anchors.sbatch`. The login node is heavily
contended (observed ~11% CPU per worker), which makes a login-node run roughly
10× slower and is antisocial besides.

## State

- `protein_ncrna/{seqio,tools,orfs,anchors,windows}.py` — steps 1–2, written and
  validated.
- `scripts/{sample_genomes,collect_anchors}.py`, `sbatch/collect_anchors.sbatch`.
- `tests/test_windows.py` — strand symmetry, contig-edge truncation, short-contig
  filter. Passing.
- `tests/test_blind.py` — config is protein-only, Rfam models rejected, no
  discovery module touches `unblind/`. Passing.
- No full-scale run submitted yet; steps 3–9 not started.
