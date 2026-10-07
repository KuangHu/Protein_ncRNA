# unblind/ — quarantine

Known RNA coordinates, Rfam models and published annotations for the three
ground-truth systems (IS110 bridge RNA, retron msr/msd, Cas9 CRISPR array and
tracrRNA) belong **here and nowhere else**.

Nothing under `protein_ncrna/` may import from this directory or read a file in
it. `tests/test_blind.py` fails the build if it does.

The unblinding comparison — `blind prediction ∩ known RNA coordinate`, scored on
overlap, rank and clade coverage — is a separate script that reads frozen
prediction output. Predictions are frozen *before* anything here is opened. If a
prediction set is regenerated after looking at this directory, it is no longer a
blind result and must not be reported as one.

Empty by design until the discovery side produces something to unblind against.
