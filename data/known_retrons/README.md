# data/known_retrons — the unblinding benchmark

**This directory is not part of blind discovery.** Nothing under
`protein_ncrna/`, and nothing that produces `selected_clades.tsv`, may read it.
It exists to answer one question *after* blind predictions are frozen: does the
pipeline recover retrons it was never told about?

See `../../unblind/README.md` for the same rule applied to the other controls.

## Sources

### RetronDB (Shipman Lab, Retron Census)

`raw/retrondb.csv` — downloaded 2026-10-06 from

```
https://raw.githubusercontent.com/Shipman-Lab/Retron-Census/main/retrondb_for_RTDNA_census.csv
```

1,928 rows. Columns: `_id`, `node`, `rt/clade`, `retron (sub)`,
`msr/msd familiy`, `cluster/domain`, `accesion`, `retron name`, `taxon code`,
`species/strain`, taxonomy, `ncRNA`.

Contents as measured:

| property | value |
|---|---|
| rows | 1,928 |
| rows carrying an ncRNA sequence | **180** |
| ncRNA length | 88 / 176 / 444 nt (min / median / max) |
| unique species or strains | 1,572 |
| RT clades | 11 (clade 8 largest at 427 rows) |
| named retrons | 17 |

### Rfam retron covariance models (optional, not yet obtained)

Only to be used under `out/known_retron_benchmark/`, and only once Infernal is
installed. Hits from an Rfam or custom msr/msd CM must never be fed back into
`selected_clades.tsv`.

## Matching is harder than it looks — read before trusting a recall number

RetronDB identifiers are **not** assembly accessions:

| accession style | rows |
|---|---|
| PATRIC `fig\|670897.3.peg.2382` | 1,577 |
| RefSeq/GenBank protein (`WP_…`, `EIJ…`) | 275 |
| unrecognised | 76 |

This pipeline calls ORFs de novo with prodigal and never assigns protein
accessions, and `genomes_db` is keyed by GCA assembly accession. **There is
therefore no accession join to make.** `annotate_known_retrons.py` still
attempts one, and it is expected to recover close to nothing; the number is
reported so the limitation stays visible rather than being assumed away.

The real route is **sequence**: blastn of the 180 known ncRNAs against extracted
RT windows. Consequences to state whenever a recall figure is quoted:

- Recall is measured against **180 ncRNAs, not 1,928 retrons**. A miss may mean
  the pipeline failed, or simply that RetronDB holds no sequence for that row.
- RetronDB is enterobacteria-heavy and so is `genomes_db`; agreement between
  them is not evidence of generality.
- A blastn hit shows a known retron ncRNA lies inside an RT window. It does not
  show the blind pipeline would have *called* that block, which is the question
  step 4 answers.

## Layout

```
data/known_retrons/
├── README.md             this file
├── raw/retrondb.csv      downloaded, never edited
└── known_retrons.tsv     normalized, written by scripts/import_retrondb.py
```

`raw/` holds untrusted downloaded data: parse it with `python -I`, never run an
interpreter from inside it.
