# configs/known_retron — benchmark configs

Unlike `configs/blind/`, files here **may** name msr, msd, retron RNA families,
Rfam accessions and covariance models. That is the point: this branch exists to
unblind predictions after they are frozen.

Nothing under `protein_ncrna/`, and none of the scripts that produce
`selected_clades.tsv`, may read this directory. `tests/test_leakage.py` fails
the build if that changes.

If Infernal is installed later, `cmsearch` of Rfam or custom msr/msd models runs
only under `out/known_retron_benchmark/`, and its hits are never fed back into
blind selection.
