# RNAMunin on the insertion's non-coding regions

Question: does RNAMunin (arXiv:2507.11950, `github.com/tnn111/RNAMunin`) rank
IS110, group II introns and IS200/605 above IS1, IS3 and IS4/5/10?

Answer: essentially no, except for group II introns, which are Rfam RF00029
and therefore inside its training set.

## What was scored

All 31,230 non-coding candidate and cluster-member sequences stage 1 emitted
from inside the insertions (`expand/nc_insertion/candidates/*.fasta`), median
101 bp. Published weights, both strands, 101 nt window, as upstream. The
TensorFlow import hangs on lr6, so the five layers were rebuilt in torch from
the `.keras` archive; an independent numpy forward pass agrees to 2e-6.

Labels from Pfam domains, one family per bag: 46 RNA-positive bags (IS110 21,
IS200/605 15, group II 10) against 147 RNA-negative (IS1 9, IS3 69,
IS4/5/10 69). AUC with a 20,000-shuffle label permutation.

## Result

| score | AUC | p |
|---|---|---|
| candidate length (stage 1) | **0.793** | <0.0001 |
| Minerva paired fraction | 0.726 | 0.006 |
| RNAMunin best 50 nt window | 0.648 | 0.001 |
| RNAMunin mean per base | 0.564 | 0.097 |

Median best-50 score by family:

```
GroupII    0.532     IS4/5/10  0.095
IS110      0.128     IS3       0.034
IS200/605  0.118     IS1       0.021
```

Group II carries the signal. Drop it and the AUC falls to 0.608; IS110 and
IS200/605 alone give the same 0.608. IS4/5/10, an RNA-negative family, scores
closer to IS110 than IS110 scores to group II.

Not a length artefact: Spearman(score, length) = 0.07. Partly a composition
artefact: Spearman(score, GC) = -0.30, and GC alone gives AUC 0.579, which is
the whole of the mean-score AUC of 0.564. Within length quartiles the AUC is
0.24 / 0.66 / 0.59 / 0.49 -- no trend, and the longest quartile, where the
padding matters least, is at chance.

Adding it to what already works makes things worse: length alone 0.793 drops
to 0.755 when multiplied by RNAMunin, and length x Minerva 0.817 drops to
0.810 when RNAMunin joins.

## Reading

RNAMunin was trained on Rfam. It recognises the family it has seen and does
not generalise to the two whose RNAs were described after Rfam closed on
them. That is a useful negative: it says the IS110 bridge RNA and the
IS200/605 omega RNA do not look like a generic Rfam ncRNA at the sequence
level, which is why a blind caller is needed for them at all.

Caveat on fit: upstream scans whole contigs and these are short fragments
padded circularly, so this is slightly off-label use. The length-quartile
breakdown is the control for that and it does not rescue the result.
