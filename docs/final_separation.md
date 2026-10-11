# RNA-guided against every other IS family: where the score stands

IS110 (bridge RNA), IS200/605 (omega RNA) and group II introns -- all three
use an RNA to find their target -- against every bag in the corpus that looks
like a transposase and has no known non-coding RNA.

## Three axes

| axis | what it measures | why it should work |
|---|---|---|
| `cand_len` | how much conserved non-coding the element carries | a guide module takes sequence |
| `minerva` | mean paired fraction over the element's non-coding regions | a guide RNA folds |
| `strand` | how consistently the candidate sits on one strand of the transposase | a transcript has a strand; a conserved stretch of DNA need not |

The third only became computable after `orient_candidates.py`, because the
insertions are stored in lexicographic orientation and the strand relative to
the protein was a coin flip until then.

Alone, over every other IS family:

```
cand_len   AUC 0.789      minerva   AUC 0.752      strand   AUC 0.703
```

Strand consistency on its own: RNA-guided median **1.000**, RNA-free 0.907
(bags with at least 3 placed window hits, n+ 37, n- 138, p = 0.0002).

## Combined

Same bags in every row, so the gain is the feature and not the subset.

| | n | AUC | IS110 | IS200/605 | GroupII |
|---|---|---|---|---|---|
| len x minerva | 175 | 0.869 | 0.838 | 0.857 | 0.937 |
| **len x minerva x strand** | 175 | **0.910** | 0.904 | 0.871 | 0.954 |
| len x minerva, IR arm < 12 | 158 | 0.888 | 0.865 | 0.874 | 0.951 |
| **len x minerva x strand, IR arm < 12** | 158 | **0.924** | 0.915 | 0.884 | 0.978 |

Adding strand is worth +0.041 raw and +0.036 after the inverted-repeat
filter, on identical bag sets both times.

## Where it came from

```
0.655   Minerva alone, aggregated by max over regions   (the aggregator was wrong)
0.819   Minerva alone, aggregated by mean
0.864   x candidate length, against every other IS family
0.895   with the inverted-repeat bags set aside
0.924   x strand consistency
```

All three families now sit above 0.88 and none of them carries the result:
drop group II and the other two are still 0.915 and 0.884.

## What is still unresolved

ISL3 remains the worst false positive and none of the three axes explains it
-- it is not inverted repeats, not unannotated protein, and it has the
largest non-coding fraction of any family measured. It is either a blind spot
or an undescribed RNA, and this benchmark cannot tell which.

RNAMunin stays out: it subtracts at every level tested.

## Committed output

`final_separation_bench.py`, output in `results/benchmarks/`. The combined
table above reproduces exactly. The single-axis line differs, because the
script measures each axis on the 175 bags that carry **all three** rather
than on every bag carrying that one axis:

```
            this script (n=175)     as first reported
cand_len         0.826                   0.789
minerva          0.754                   0.752
strand           0.685                   0.703
```

The common-bag version is the one to quote: it is the set the combinations
are scored on, so the axes and their combinations are then comparable.
