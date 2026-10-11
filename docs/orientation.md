# Is the insertion standardized to the protein's direction?

No. Every one of the 78,120 sites in `insertions_sites.jsonl` carries
`orient_source: canonical_lexicographic_min` -- the insertion was flipped to
whichever of the sequence and its reverse complement sorts first. That is
deterministic, and it is a coin flip against the transposase:

```
orient:  rc 39,100   fwd 39,020
```

1,715 of 6,139 bags hold sites in both orientations. Blasting each bag's
candidate against its own cluster members, **5,230 homologous member regions
match only on the minus strand** -- 31% of all homologous members are stored
reverse-complemented relative to their own candidate.

## What this does not break

Checked rather than assumed, because the obvious worry is that the clusters
were split in half:

- **mmseqs clusters a sequence with its reverse complement.** Direct test: a
  random 300 bp sequence, its reverse complement and an unrelated sequence,
  `easy-cluster` with the pipeline's settings puts `fwd` and `rev` in one
  cluster. So recurrence, distinct-site counts and `candidate_len` are
  unaffected.
- **Windows are already turned to the anchor.** `windows.extract_window` does
  `orient_to_anchor_strand`, so every Minerva region is in the transposase's
  frame.
- **The side call is read from the windows, not from the insertion record.**
  `element_geometry.py` blasts the non-coding regions into the anchor-oriented
  windows and assigns upstream/downstream in the frame where the tnp occupies
  `[0, anchor_len)`. So `[nc][tnp]` against `[tnp][nc]` is correct: 198 bags
  both sides, 48 `[nc][tnp]`, 20 `[tnp][nc]`, 124 unplaced.
- **Every score reported is strand-symmetric.** Minerva folds both strands,
  RNAMunin scores both strands, the inverted-repeat arm is symmetric by
  definition and the ORF z-score scans six frames.

So no AUC in any benchmark moves because of this.

## What it does break

The candidate sequence itself. Half the candidate FASTAs are written
3'->5' relative to their element, so the file we would hand to anyone is a
coin flip on strand. Anything downstream that cares which strand is
transcribed -- promoter, terminator, the guide portion of a bridge-RNA-like
module, a probe or a primer -- cannot be read off it.

## The fix

`scripts/orient_candidates.py` maps each candidate back into its own
anchor-oriented windows and re-emits it on the strand the window shows,
with the side it sits on, into `nc_insertion/oriented/`. It writes new files
and changes no score, so the pipeline stays frozen.

## Result

333 candidates mapped back into their windows; 309 have no windows and 41 no
hit above 90% identity.

```
flipped 160 / 333 = 48%
```

A coin flip, as predicted. Now that the frame is the transposase's, the side
each candidate sits on is readable, and it is not uniform:

| group | n | upstream | downstream | overlaps tnp |
|---|---|---|---|---|
| IS1 | 9 | **89%** | 11% | 0% |
| IS3 | 44 | **70%** | 20% | 9% |
| ISL3 | 24 | 58% | 38% | 4% |
| IS110 | 21 | 43% | 48% | 10% |
| IS200/605 | 9 | 44% | 44% | 11% |
| GroupII | 10 | 10% | **80%** | 10% |
| IS4/5/10 | 52 | 40% | 35% | **25%** |
| all | 333 | 50% | 35% | 16% |

IS1 and IS3 put their conserved non-coding almost entirely upstream of the
transposase; group II puts it downstream; IS4/5/10 has by far the most
candidates overlapping the ORF, which fits the earlier finding that those
elements carry transposase ORFs far shorter than their element.

## The Ldr check

Both type I toxin-antitoxin bags come back the same way:

```
CDS00749  377 bp  revcomp  overlaps_tnp  frac_modal_side 1.000
CDS01102  285 bp  revcomp  overlaps_tnp  frac_modal_side 1.000
```

Every placed hit puts the candidate **on top of the toxin ORF**, and both
needed flipping. That is the architecture of a type I antitoxin: RdlD is
transcribed antisense to, and overlapping, ldrD. The blind candidate landed
on the right sequence in the right relationship to the protein, and only the
orientation pass makes that visible.
