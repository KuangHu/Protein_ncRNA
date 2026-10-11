# Minerva separates these families; the max aggregator hid it

Scoring a bag by the **maximum** `minerva_paired_fraction` over its non-coding
regions gives AUC 0.655. Scoring it by the **mean** gives **0.819**, on the
same 42 RNA-positive and 137 RNA-negative bags.

The max is a draw from the tail. A bag with twenty regions gets twenty
chances at a lucky fold, so the max partly measures how many regions a bag
has; the mean asks the question that was meant -- are this element's
non-coding regions structured.

| aggregation | n+ / n- | AUC | IS110 | IS200/605 | GroupII |
|---|---|---|---|---|---|
| max over regions | 42 / 137 | 0.655 | 0.545 | 0.720 | 0.812 |
| mean over regions | 42 / 137 | **0.819** | 0.750 | 0.885 | 0.889 |
| mean, regions >= 200 bp | 41 / 101 | **0.878** | 0.839 | 0.887 | 0.946 |

## Length is a floor on folding, not a confound

Region-level AUC is 0.707 over all 3,161 labelled regions, and it rises
monotonically with how much sequence there is to fold:

```
region length   40-106   106-175   175-300   300-1503
AUC              0.653     0.614     0.728      0.815
```

Spearman(paired fraction, region length) is only 0.11, so this is not the
score tracking length. A folding model needs a fold's worth of sequence.

Within the regions of at least 300 bp, against the same long RNA-negative
regions: IS110 0.774, IS200/605 0.826, group II 0.866, and IS110 plus
IS200/605 with group II removed entirely still 0.789. No single family
carries it.

## Combined with the stage-1 candidate length

| score | n | AUC | IS110 | IS200/605 | GroupII |
|---|---|---|---|---|---|
| candidate length | 193 | 0.793 | 0.879 | 0.634 | 0.849 |
| Minerva mean | 179 | 0.819 | 0.750 | 0.885 | 0.889 |
| **length x Minerva mean** | 177 | **0.900** | 0.886 | 0.881 | 0.949 |

The complementarity is real and it is now measured on 42 RNA-positive bags
rather than the 11 of the overlap run: length carries IS110, Minerva carries
IS200/605, and the product carries both.
