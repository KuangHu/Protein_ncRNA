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

## The three negative families are not equally negative

| family | n | candidate len (med) | Minerva mean (med) | combined rank pct (med) | above a cutoff keeping 90% of RNA+ |
|---|---|---|---|---|---|
| IS1 | 9 | 86 | 0.078 | 0.33 | 0 / 9 (0%) |
| IS3 | 65 | 84 | 0.135 | 0.34 | 11 / 65 (17%) |
| IS4/5/10 | 59 | 119 | 0.157 | 0.49 | 25 / 59 (42%) |
| IS200/605 | 11 | 125 | 0.290 | 0.77 | - |
| IS110 | 21 | 250 | 0.195 | 0.75 | - |
| GroupII | 10 | 245 | 0.277 | 0.82 | - |

IS4/5/10 sits significantly above IS1 and IS3: AUC 0.702, p < 0.0001. The gap
is **length**, not structure -- 0.723 on candidate length against 0.581
(p = 0.054) on the Minerva mean. Inside IS4/5/10 the elevated bags are the
ones carrying DUF772 or DDE_Tnp_1_6, the IS5/IS1182 accessory domains:
median 132-162 bp of non-coding against 107-112 bp without them.

So the headline number depends on which negative you ask about:

```
RNA+ vs IS1 + IS3      AUC 0.949
RNA+ vs all RNA-       AUC 0.898
RNA+ vs IS4/5/10       AUC 0.835
```

One label caveat worth stating: IS10 carries RNA-OUT, a ~70 nt antisense
regulator of transposase translation, so "IS4/5/10 has no ncRNA" is not
strictly true. The data does not support RNA-OUT explaining the elevation --
RNA-OUT is short and structured while this gap is long and unstructured, and
the top IS4/5/10 outliers hold 262-356 bp of non-coding, not a 70 nt stem
loop. Those 25 bags above the cutoff are the obvious thing to look at next.
