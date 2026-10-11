#!/usr/bin/env python
"""Score sequences with RNAMunin (Lui & Nielsen 2025, arXiv:2507.11950).

The upstream `predict13` loads the model, reads one FASTA record and prints a
per-base probability for each strand. That is the right interface for one
contig and the wrong one for 31k short non-coding regions, so this does the
same arithmetic in batch: same 101 nt window, same one-hot order, same
circular padding, same forward/reverse pair. The model file and its weights
are untouched.

Per sequence it reports the summaries that could plausibly separate an RNA
from a random intergenic stretch: the mean per-base probability (length
independent by construction), the best 50 nt window, and the fraction of
bases over 0.5. Scoring only the core -- the bases whose 101 nt window lies
entirely inside the sequence -- is reported alongside so a result cannot be
blamed on the wrap-around padding.
"""
import argparse
import os
import sys

import numpy as np

SEQ_LENGTH = 101
HALF = SEQ_LENGTH // 2
MAP = {"A": 0, "T": 1, "G": 2, "C": 3, "a": 0, "t": 1, "g": 2, "c": 3}
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def read_fasta(path):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def onehot(seq):
    """N and anything else become an all-zero column, not a fake base."""
    a = np.zeros((len(seq), 4), dtype=np.float32)
    for i, c in enumerate(seq):
        j = MAP.get(c)
        if j is not None:
            a[i, j] = 1.0
    return a


def pad(seq):
    """Circular padding, as upstream does for a chromosome. A fragment
    shorter than the pad is tiled until it is long enough."""
    s = seq
    while len(s) < HALF:
        s = s + s
    return s[-HALF:] + seq + s[:HALF]


def windows(seq):
    f = onehot(pad(seq))
    n = len(seq)
    idx = np.arange(n)[:, None] + np.arange(SEQ_LENGTH)[None, :]
    return f[idx]


def summarise(p, n):
    """p: per-base probability, 5' to 3' of the given sequence."""
    out = {"mean": float(p.mean()), "max": float(p.max()),
           "frac_gt50": float((p > 0.5).mean())}
    w = min(50, n)
    if n >= w:
        c = np.concatenate([[0.0], np.cumsum(p)])
        roll = (c[w:] - c[:-w]) / w
        out["best50"] = float(roll.max())
    else:
        out["best50"] = out["mean"]
    if n > SEQ_LENGTH:
        core = p[HALF:n - HALF]
        out["core_mean"] = float(core.mean())
        out["core_n"] = len(core)
    else:
        out["core_mean"] = float("nan")
        out["core_n"] = 0
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta")
    ap.add_argument("--model", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch", type=int, default=4096)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
    import keras
    keras.config.disable_interactive_logging()
    model = keras.saving.load_model(a.model, compile=False)

    recs = list(read_fasta(a.fasta))
    if a.limit:
        recs = recs[:a.limit]
    sys.stderr.write("%d sequences\n" % len(recs))

    fields = ["mean", "best50", "max", "frac_gt50", "core_mean"]
    with open(a.out, "w") as out:
        out.write("seq_id\tbag\trole\tlen\t" +
                  "\t".join(f + s for s in ("_fwd", "_rc") for f in fields) +
                  "\tcore_n\n")
        buf, meta = [], []
        done = 0

        def flush():
            nonlocal buf, meta, done
            if not buf:
                return
            x = np.concatenate(buf, axis=0)
            p = model.predict(x, batch_size=a.batch, verbose=0).reshape(-1)
            k = 0
            for name, seq, nf in meta:
                n = len(seq)
                f = summarise(p[k:k + n], n)
                k += n
                r = summarise(p[k:k + n][::-1], n)   # rc scored 3'->5'
                k += n
                parts = name.split("|")
                bag = parts[0]
                role = parts[1] if len(parts) > 1 else "."
                out.write("%s\t%s\t%s\t%d\t%s\t%d\n" % (
                    name, bag, role, n,
                    "\t".join("%.6f" % f[q] for q in fields) + "\t" +
                    "\t".join("%.6f" % r[q] for q in fields), f["core_n"]))
            done += len(meta)
            sys.stderr.write("%d\n" % done)
            sys.stderr.flush()
            buf, meta = [], []

        pend = 0
        for name, seq in recs:
            w = windows(seq)
            rc = windows(seq.translate(COMP)[::-1])
            buf.append(w)
            buf.append(rc)
            meta.append((name, seq, len(seq)))
            pend += 2 * len(seq)
            if pend >= 200000:
                flush()
                pend = 0
        flush()


if __name__ == "__main__":
    main()
