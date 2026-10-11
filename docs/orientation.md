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
