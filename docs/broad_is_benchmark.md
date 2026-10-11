# The three RNA carriers against every other IS family

The six-family benchmark asked a narrow question. This widens the negative to
every bag in the corpus that looks like a transposase and has no known
non-coding RNA: ISL3, IS256, IS66, Mu, Tn3 and IS630 on top of IS1, IS3 and
IS4/5/10. 42 RNA-positive bags against **176** negatives.

Widening a negative class is where a benchmark gets quietly corrupted, so a
bag is excluded rather than counted as a negative if it carries a reverse
transcriptase domain (retrons, DGRs, group II), a Cas domain, TnsB (Tn7-like
elements are guided by a crRNA), HNH (the HEARO RNAs) or a phage integrase.
Three bags were dropped that way, 14 more for hitting two families at once.

## Result

| negative set | n | n- | AUC | p |
|---|---|---|---|---|
| IS1 + IS3 + IS4/5/10 | 163 | 121 | 0.890 | <0.0001 |
| **all other IS** | 218 | 176 | **0.864** | <0.0001 |
| all other IS, IR arm < 12 | 191 | 153 | **0.895** | <0.0001 |

Per family, against the full negative set: IS110 0.841, IS200/605 0.848,
group II 0.930. With the inverted-repeat bags removed: 0.879, 0.878, 0.957.

So the separation survives the harder negative. It costs about 0.03.

## Which negatives misbehave

Median combined-rank percentile, and how many clear a cutoff that keeps 90%
of the RNA-positive bags:

| family | n | med rank | above cutoff | med Minerva |
|---|---|---|---|---|
| IS1 | 9 | 0.29 | 0 / 9 | 0.078 |
| IS3 | 52 | 0.32 | 10 / 52 | 0.152 |
| IS256 | 10 | 0.36 | 0 / 10 | 0.158 |
| IS66 | 9 | 0.36 | 3 / 9 | 0.138 |
| IS630 | 3 | 0.27 | 1 / 3 | 0.162 |
| Mu | 3 | 0.51 | 0 / 3 | 0.126 |
| IS4/5/10 | 60 | 0.49 | 25 / 60 | 0.156 |
| Tn3 | 5 | 0.62 | 3 / 5 | 0.175 |
| **ISL3** | **25** | **0.66** | **21 / 25** | **0.260** |

The families added here are mostly well behaved. IS256, IS66 and IS630 sit
with IS1 and IS3. ISL3 does not.

## ISL3 is the worst false positive in the benchmark

ISL3 scores **higher than IS4/5/10 and its median Minerva score, 0.260, is
above IS110's 0.195** -- third of any family here, behind only IS200/605
(0.290) and group II (0.276). 21 of its 25 bags clear the cutoff.

Neither explanation that worked for IS4/5/10 works here:

| | ISL3 | IS4/5/10 | IS110 |
|---|---|---|---|
| median IR arm | 8 | 8 | 9 |
| frac IR arm >= 12 | 0.12 | 0.24 | 0.05 |
| median ORF z vs shuffle | 0.64 | 0.19 | 0.35 |
| median element / ORF (bp) | 1508 / 852 | 1348 / 924 | 1432 / 1005 |
| non-coding fraction | **0.43** | 0.31 | 0.30 |

Not inverted repeats -- the IR filter moves ISL3 from 0.66 to 0.67. Not
unannotated protein -- the ORF z-score is modest. What ISL3 does have is the
largest non-coding fraction of any family measured, 43% of the element, and
that non-coding folds.

It is not one recurring artefact either: 18 of 31 ISL3 bags share candidate
sequence with another bag, but across 4 separate components, and 13 bags are
singletons.

## How to read this

Two possibilities, and this benchmark cannot choose between them. Either the
score has a blind spot on ISL3 elements, or ISL3 carries a structured
non-coding RNA that is not in the literature. The second is the thing this
pipeline was built to find, so it is a lead, not a conclusion -- and the
honest statement today is that ISL3 is an unexplained positive, measured
blind, with no inverted repeat and no hidden ORF behind it.

Nothing here feeds the caller. `scripts/broad_is_bench.py` holds the family
map and the exclusion list; `data/broad_is_labels.tsv` the labels.
