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

| score | n+ / n- | AUC | p |
|---|---|---|---|
| candidate length (stage 1) | 46 / 147 | **0.793** | <0.0001 |
| Minerva paired fraction, mean over a bag's regions | 42 / 137 | 0.819 | <0.0001 |
| Minerva paired fraction, max over a bag's regions | 42 / 137 | 0.655 | 0.001 |
| RNAMunin best 50 nt window | 46 / 145 | 0.648 | 0.001 |
| RNAMunin mean per base | 46 / 145 | 0.564 | 0.097 |

Length x Minerva, on the 177 bags carrying both, gives **AUC 0.900** with
IS110 0.886, IS200/605 0.881 and group II 0.949 -- the first score that is
not carried by one family. RNAMunin added as a third term drops it to 0.844.

The max over a bag's regions was the wrong aggregator and cost 0.16 of AUC:
it scores a bag on its single luckiest fold, which is a draw from the noise,
while the mean asks whether the bag's non-coding regions are structured as a
population. Restricting to regions of at least 200 bp -- enough sequence for
a fold to exist -- raises Minerva alone to 0.878.

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
to 0.755 when multiplied by RNAMunin, and length x Minerva 0.900 drops to
0.844 when RNAMunin joins.

## Reading

RNAMunin was trained on Rfam. It recognises the family it has seen and does
not generalise to the two whose RNAs were described after Rfam closed on
them. That is a useful negative: it says the IS110 bridge RNA and the
IS200/605 omega RNA do not look like a generic Rfam ncRNA at the sequence
level, which is why a blind caller is needed for them at all.

Caveat on fit: upstream scans whole contigs and these are short fragments
padded circularly, so this is slightly off-label use. The length-quartile
breakdown is the control for that and it does not rescue the result.

## Re-tested with the treatment Minerva got

Minerva only worked once it was aggregated by the mean over a bag's regions
and restricted to regions long enough to fold. RNAMunin was first scored
without either, so it was re-run across 20 settings: mean and max, best-50
and per-base mean, at minimum region lengths of 0, 100, 150, 200 and 300 bp.

Best cell: minimum length 0, best-50, mean over regions, **AUC 0.656**. The
length cutoff makes RNAMunin *worse* (0.656 -> 0.565 at 150 bp -> 0.546 at
200 bp) where it made Minerva better (0.819 -> 0.878). That contrast is
itself diagnostic: Minerva improves when given enough sequence to fold
because it is reading a fold; RNAMunin does not, because it is not.

### One correction

Saying it "only sees group II introns" was too strong. Under the better
aggregation it carries a weak IS110 signal as well:

| family | RNAMunin | Minerva |
|---|---|---|
| IS110 | 0.660 | 0.744 |
| IS200/605 | 0.528 | 0.882 |
| GroupII | 0.798 | 0.886 |

IS200/605 is at chance. Group II is still the strongest, which is still the
family in its training set. And it is dominated by Minerva on all three.

### It is also reading the inverted repeats

Spearman(IR arm, RNAMunin best-50) = 0.254, against 0.361 for Minerva. So
the hairpin confound is not specific to the folding model; a sequence CNN
trained on Rfam picks up transposon ends too. Removing the long-IR bags
moves RNAMunin only from 0.658 to 0.669.

### It never helps the combination

| bags | length x Minerva | + RNAMunin |
|---|---|---|
| all 175 | 0.898 | 0.855 |
| IR arm < 12 (154) | **0.934** | 0.893 |

Tested at both the raw and the IR-cleaned level, with the best aggregation
RNAMunin has. It subtracts in both. It stays out.
