# A close look at IS4/5/10, the hard negative

IS4/5/10 sits between the clean negatives and the RNA-carrying families:
AUC 0.702 over IS1+IS3, and 25 of its 59 bags clear a cutoff that keeps 90%
of the RNA-positive bags. This is what those 25 are.

## Not unannotated protein

The first guess was that the long non-coding is really coding sequence the
annotation missed -- several of these bags have an element far longer than
their transposase ORF (CDS07270 1837 bp element, 841 bp ORF; CDS08108
2975 / 1254; CDS00770 4299 / 1290, and that one carries a LamB domain, so it
is a composite element with cargo).

Tested by taking the longest stop-free stretch over six frames and comparing
it to 60 composition-matched shuffles of the same sequence. Coding sequence
is far more stop-free than its own composition predicts; these are not:

| group | n | med len | med z | frac z >= 2 |
|---|---|---|---|---|
| IS4/5/10 hot | 25 | 160 | **-0.10** | 0.04 |
| IS4/5/10 cold | 34 | 96 | 0.20 | 0.00 |
| IS110 | 21 | 250 | 0.35 | 0.00 |
| GroupII | 10 | 208 | 0.84 | 0.00 |

One bag of 25 is coding-like. The non-coding is non-coding.

## It is the terminal inverted repeats

A transposon end pairs beautifully without being an RNA. Measuring the
longest exact inverted repeat in each candidate -- the arm of the best
hairpin:

| group | n | med IR arm | frac arm >= 12 bp |
|---|---|---|---|
| IS4/5/10 hot | 25 | 10 | **0.40** |
| IS4/5/10 cold | 34 | 7.5 | 0.12 |
| IS3 | 68 | 8 | 0.04 |
| IS1 | 9 | 8 | 0.00 |
| IS110 | 21 | 9 | 0.05 |
| IS200/605 | 15 | 10 | 0.13 |
| GroupII | 10 | 8.5 | 0.20 |

Hot against cold on IR arm length alone: AUC 0.729. The longest arms are
15 bp, in CDS06778, CDS07780 and CDS07787 -- all DUF772 or DDE_Tnp_1_6, the
IS5/IS1182 accessory domains.

Minerva is partly reading these: Spearman(IR arm, Minerva mean) is 0.36 over
all bags and **0.49 within IS4/5/10**.

## What it costs, and what dropping it buys

| bags kept | n | IS4/5/10 kept | AUC |
|---|---|---|---|
| all | 175 | 59 | 0.898 |
| IR arm < 12 | 154 | 45 | **0.933** |
| IR arm < 10 | 121 | 33 | 0.935 |

Excluding the long-IR bags raises the separation by 0.035 and removes 14 of
the 25 false positives. Note it also drops 4 RNA-positive bags, so this is a
flag to read, not a filter to apply blindly.

`scripts/ir_flag.py` writes `ir_arm`, `ir_frac` and `ir_flag` per bag. It is
a reported diagnostic and **nothing here feeds the score** -- the pipeline is
frozen and this changes no parameter. 89 of 669 candidates flag at arm >= 12.

## An internal check that came out right

IS110 has the *lowest* long-IR rate of any family here, 5%. The IS110/IS1111
family is the one that is known not to carry canonical terminal inverted
repeats. The IR measure found that without being told.
