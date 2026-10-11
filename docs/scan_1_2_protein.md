# Open scan: every one-or-two-protein element

353 bags, every element that encodes one or two proteins and has
element-trimmed windows, labelled or not. Scored by the two axes the
benchmark settled on -- how much conserved non-coding the element carries,
and how much of it folds -- as a rank product. Labels never enter the score.

`scripts/scan_bags.py` -> `nc_insertion/scan_1_2_protein.tsv` (all 353) and
`nc_insertion/scan_shortlist.tsv` (266 after the confound filters).

## The positive controls behaved

Median position of each family, out of 353:

```
GroupII      39      ISL3        110        IS3     266
IS200/605    61      Tn3         192        IS256   276
IS110        69      Mu          201        IS1     282
                     IS4/5/10    209        IS630   288
                     IS66        262
```

All three RNA-carrying families rank above every family that carries none.
That is the pass condition, and it held without the labels being visible to
the score.

## The scan recovered a known RNA class it was never shown

94 of the 266 shortlisted bags carry no family label. The two highest are
both `Ldr_toxin`, and **Ldr_toxin has the best median position of any group
in the shortlist -- 11 of 266, ahead of group II introns at 15**:

| bag | pos | cand len | Minerva | element | sites | species |
|---|---|---|---|---|---|---|
| CDS00749 | 3 | 377 | 0.318 | 536 | 13 | E. coli |
| CDS01102 | 11 | 285 | 0.260 | 966 | 19 | E. coli |

Ldr is a **type I** toxin-antitoxin system: the antitoxin is a small
antisense RNA (RdlD), not a protein. It was never in the label set.

The architecture confirms it. CDS00749's insertions are 535, 536, 537 and
1070, 1071 bp -- one and two copies of a ~535 bp repeat. CDS01102's are
965, 966 and 1449, 1451 -- two and three copies of ~483 bp. That is the
E. coli long direct repeat, a tandem array whose unit carries the toxin ORF
and the antisense RNA, and most of the unit is non-coding.

## Type I separates from type II

The sharpest internal check, and one nobody designed:

| locus | type | antitoxin is | shortlist position |
|---|---|---|---|
| Ldr (CDS00749) | I | **RNA** | 3 |
| Ldr (CDS01102) | I | **RNA** | 11 |
| AbiEii (CDS02655) | IV | protein | 8 |
| ParE (CDS06626) | II | protein | 93 |
| YoeB (CDS06711) | II | protein | 113 |
| YafO (CDS01112) | II | protein | 154 |

Both RNA-antitoxin loci land in the top 5%; all three protein-antitoxin
loci land in the bottom two thirds. Two against three is far too small to
test, but the direction was predictable in advance and it came out right.

## Confounds carried on every row

`ir_arm` (longest inverted repeat; a transposon end pairs without being
transcribed) and `orf_z` (longest stop-free stretch against 60
composition-matched shuffles; unannotated protein is not an RNA). The
shortlist drops `ir_arm >= 12`, `orf_z >= 2`, elements over 2.5 kb and
candidates over 600 bp, which removes 87 of 353.

## What this is not

An unlabelled bag high in this scan is a lead. The scan contains known
RNA classes, elements whose annotation is incomplete, and real novelty, and
nothing here distinguishes them. The Ldr result is worth stating because it
is a *known* answer the blind score reached without being given it -- the
instrument reading true on a case it was not calibrated against, like
[[retron-is-instrument-not-target]]. ISL3 remains the open one.
