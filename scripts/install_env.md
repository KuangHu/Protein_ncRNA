# Environment

## Steps 1–2 and the census need only what is already installed

```
prodigal      /global/home/users/kh36969/.conda/envs/opfi/bin/prodigal
hmmsearch     (same env)
hmmfetch      (same env)
mmseqs        (same env)   # divergence ladder in census_report.py
python 3.11   (same env)
```

Pfam-A: `/global/scratch/users/kh36969/pfam_db/Pfam-A.hmm` (pressed; `hmmfetch`
writes its index as `Pfam-A.hmm.h3m.ssi`).

Binaries are resolved from the directory of the running interpreter
(`protein_ncrna/tools.py`), so an absolute `$PY` path works unactivated.

**Do not install the RNA tools before the census.** Nothing in steps 1–2 or the
census touches them.

## R-scape and Infernal — required later, deliberately deferred

They are needed only at step 6, conserved noncoding block → covariation →
covariance model. Installing them first risks the embarrassing outcome of a
working R-scape with no clade that has enough informative substitutions to feed
it. The order is:

```
anchor census
  → pick clades with real diversity (census_report.py: covariation_feasible)
    → install RNA tools
      → step 6
```

When that point is reached, build a separate env rather than perturbing `opfi`,
which other projects depend on:

```bash
conda create -n rna_cov -c conda-forge -c bioconda infernal rscape
conda activate rna_cov && cmbuild -h && R-scape -h
```

Already present in `opfi` for step 6's earlier stages: `mafft`, `RNAalifold`,
`RNAfold`, `esl-alistat`, `blastn`.

## genomes_db is 74% empty stubs

Do not scan `INDEX.tsv` directly. Generate the manifest once:

```bash
$PY scripts/genomes_db_report.py --out out --verify-gzip   # SLURM, reads every file
```

and treat `out/real_genomes.tsv` as the shared entry point for every pipeline
that reads this mirror. A size filter alone (`--min-bytes`, the default) is fast
and catches the 20-byte stubs; `--verify-gzip` additionally catches truncated
downloads that pass on size.
