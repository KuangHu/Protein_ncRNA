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



# --- the model -------------------------------------------------------------
#
# Convolution13.keras is a plain Sequential stack and TensorFlow will not
# import on this cluster's compute nodes, so the weights are read straight out
# of the .keras archive and run under torch. Keras Conv1D stores a kernel as
# (width, in, out) and Flatten reads channels-last in C order; both are
# reproduced here, and --selftest re-runs the first batch through an
# independent numpy implementation so a layout mistake cannot pass silently.

LAYERS = ["conv1d", "conv1d_1", "conv1d_2", "dense", "dense_1"]


def load_weights(path):
    import io
    import zipfile
    import h5py
    with zipfile.ZipFile(path) as z:
        blob = z.read("model.weights.h5")
    f = h5py.File(io.BytesIO(blob), "r")
    w = {}
    for name in LAYERS:
        g = f["layers"][name]["vars"]
        w[name] = (np.array(g["0"]), np.array(g["1"]))
    return w


def numpy_forward(x, w):
    def conv(a, k, b):
        n, L, cin = a.shape
        kw, _, cout = k.shape
        out = np.zeros((n, L - kw + 1, cout), dtype=np.float32)
        flat = k.reshape(kw * cin, cout)
        for i in range(L - kw + 1):
            out[:, i, :] = a[:, i:i + kw, :].reshape(n, -1) @ flat
        return np.maximum(out + b, 0.0)

    def pool(a, s):
        n, L, c = a.shape
        L -= L % s
        return a[:, :L, :].reshape(n, L // s, s, c).max(axis=2)

    h = conv(x, *w["conv1d"])
    h = pool(h, 4)
    h = conv(h, *w["conv1d_1"])
    h = pool(h, 3)
    h = conv(h, *w["conv1d_2"])
    h = h.reshape(h.shape[0], -1)
    h = np.maximum(h @ w["dense"][0] + w["dense"][1], 0.0)
    z = h @ w["dense_1"][0] + w["dense_1"][1]
    return (1.0 / (1.0 + np.exp(-z))).reshape(-1)


def load_backend(path, backend, threads):
    if backend == "keras":
        os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
        import keras
        keras.config.disable_interactive_logging()
        model = keras.saving.load_model(path, compile=False)

        def run(x, batch):
            return model.predict(x, batch_size=batch, verbose=0).reshape(-1)
        return run

    import torch
    import torch.nn as nn
    if threads:
        torch.set_num_threads(threads)
    w = load_weights(path)

    def conv(name):
        k, b = w[name]
        kw, cin, cout = k.shape
        c = nn.Conv1d(cin, cout, kw)
        c.weight.data = torch.from_numpy(k.transpose(2, 1, 0).copy())
        c.bias.data = torch.from_numpy(b.copy())
        return c

    def dense(name):
        k, b = w[name]
        d = nn.Linear(k.shape[0], k.shape[1])
        d.weight.data = torch.from_numpy(k.T.copy())
        d.bias.data = torch.from_numpy(b.copy())
        return d

    c0, c1, c2 = conv("conv1d"), conv("conv1d_1"), conv("conv1d_2")
    d0, d1 = dense("dense"), dense("dense_1")
    net = nn.Sequential(c0, nn.ReLU(), nn.MaxPool1d(4),
                        c1, nn.ReLU(), nn.MaxPool1d(3),
                        c2, nn.ReLU())
    net.eval()

    def run(x, batch):
        out = []
        with torch.no_grad():
            for i in range(0, len(x), batch):
                t = torch.from_numpy(x[i:i + batch]).permute(0, 2, 1)
                h = net(t).permute(0, 2, 1).reshape(t.shape[0], -1)
                h = torch.relu(d0(h))
                out.append(torch.sigmoid(d1(h)).reshape(-1).numpy())
        return np.concatenate(out) if out else np.zeros(0, dtype=np.float32)

    run.weights = w
    return run


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta")
    ap.add_argument("--model", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch", type=int, default=4096)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--backend", choices=("torch", "keras"), default="torch")
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--selftest", action="store_true",
                    help="run the numpy reference on the first batch too")
    a = ap.parse_args()

    predict = load_backend(a.model, a.backend, a.threads)

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
            p = predict(x, a.batch)
            if a.selftest and done == 0 and hasattr(predict, "weights"):
                ref = numpy_forward(x[:256], predict.weights)
                d = float(np.abs(ref - p[:256]).max())
                sys.stderr.write("selftest max |torch-numpy| = %.3g\n" % d)
                if d > 1e-4:
                    raise SystemExit("backend mismatch: %g" % d)
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
