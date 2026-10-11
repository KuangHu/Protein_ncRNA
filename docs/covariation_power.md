# The covariation model: right test, no power on this corpus

The third opinion the repo already carries is `scripts/fold_blocks.py` --
MAFFT, RNAalifold, R-scape. It is the one test an inverted repeat should
fail: an IS end is conserved, but it has no reason to *covary*, because the
constraint on it is DNA binding, not a fold. After finding that Minerva's
IS4/5/10 false positives are terminal inverted repeats, covariation is the
obvious next discriminator.

It cannot be run meaningfully here, and the measurement says why.

## Within-bag identity

Per bag, the candidate against every member of its cluster, blastn, identity
over the shorter sequence:

| family | n | median mean-identity | lowest bag | frac of bags under 95% |
|---|---|---|---|---|
| IS110 | 20 | 99.6 | 85.2 | 0.15 |
| IS200/605 | 12 | 98.9 | 94.2 | 0.17 |
| GroupII | 7 | 99.7 | 76.4 | **0.43** |
| IS1 | 8 | 97.7 | 93.8 | 0.12 |
| IS3 | 54 | 99.7 | 83.4 | 0.15 |
| IS4/5/10 | 56 | 99.6 | 94.2 | 0.02 |

**Median 99.5% across all 157 labelled bags. 18 bags (11%) fall under 95%,
10 under 90%.** R-scape needs substitutions to test; with none it returns
"no significant covariation" for want of evidence rather than for want of
structure. `fold_blocks.py` already labels that case `underpowered` instead
of negative -- it would label roughly 89% of these bags that way.

This is the same fact that makes within-bag recurrence useless at AUC 0.608:
a bag is a pile of near-identical copies of one element. One corpus
property, two dead metrics.

## Pooling across bags does not rescue it

All-against-all blastn of the 191 labelled candidates gives 202 cross-bag
hits, which form 8 components of 3 or more bags covering 38 of 191 bags.
Only **one** of the 8 has a median cross-bag identity below 95%.

```
comp size families            med identity
   1    9  IS4/5/10:9           100.0
   2    7  IS1:7                 94.2   <- the only diverged one
   3    5  IS4/5/10:5            96.9
   4    4  IS200/605:3, IS3:1    98.2
   5    4  IS4/5/10:4            95.8
```

Component 3 is CDS06778, CDS07270, CDS07341, CDS07780 and CDS07787 -- four
of the five are the long-IR bags from `is4_5_10_closeup.md`. So those 25
IS4/5/10 false positives are not 25 independent observations; they collapse
into a handful of recurring subfamilies, and the IS4/5/10 false-positive
count should be read accordingly.

Group II is the only family with real depth, 43% of its bags under 95%,
which is what a element mobile across diverged hosts should look like.

## What would change the answer

More divergence, not another model. The corpus is nine species with most
bags drawn from one of them; covariation needs homologous non-coding regions
from genera far enough apart to have substituted. Until then the covariation
verdict on almost every bag would be `underpowered`, and reporting that as
"no structure" would be the silent failure the script was written to avoid.

A qualitative readout is still available on the 18 powered bags -- 8 RNA+
against 10 RNA- -- which is too small for an AUC but large enough to see
whether the direction is right. That needs R-scape and MAFFT installed; the
`rnafold` env named in `sbatch/install_rna_tools.sbatch` is not on this
system, and `vrna` carries RNAalifold but neither R-scape nor MAFFT.

## Caveat on the identity measure

The `memberN` records in a candidate FASTA are every non-coding region of
each site in the cluster, not only the clustered region -- `[nc][tnp][nc]`
contributes both sides. So about a third of member records align to the
candidate, and identity is measured over the homologous subset by
construction. That is the right denominator for this question, but it means
the numbers describe the cluster, not the bag's whole non-coding content.
