#!/usr/bin/env python
"""Put each non-coding candidate into the transposase's reading frame.

`insertions_sites.jsonl` stores every insertion in `canonical_lexicographic_min`
orientation: whichever of the sequence and its reverse complement sorts first.
That is deterministic but it is a coin flip against the protein -- 39,100 of
78,120 sites are stored reverse-complemented -- so a candidate's stored strand
says nothing about which strand is transcribed.

Nothing measured so far depends on it. mmseqs clusters a sequence with its
reverse complement, the windows are already turned to the anchor's strand, and
Minerva, RNAMunin, the inverted-repeat arm and the ORF z-score are all
strand-symmetric. What depends on it is handing the sequence to anyone: half
the candidates are written 3'->5' relative to the element.

So the candidate is mapped back into its own anchor-oriented windows, where the
transposase occupies [0, anchor_len), and re-emitted on the strand the window
shows, with the side it sits on. This writes new files and changes no score.
"""
import argparse
import glob
import json
import os
import subprocess
import sys
import tempfile
from collections import Counter

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(COMP)[::-1]


def read_candidate(path):
    keep, buf = False, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                keep = "|candidate|" in line
            elif keep:
                buf.append(line)
    return "".join(buf)


def window_meta(jsonl):
    m = {}
    with open(jsonl) as fh:
        for line in fh:
            d = json.loads(line)
            m[d["anchor_id"]] = (d.get("anchor_offset", 0), d.get("anchor_len", 0))
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--windows", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--blastn", default="/global/home/users/kh36969/.conda/"
                                        "envs/opfi/bin/blastn")
    ap.add_argument("--makeblastdb",
                    default="/global/home/users/kh36969/.conda/envs/opfi/bin/"
                            "makeblastdb")
    ap.add_argument("--min-identity", type=float, default=90.0)
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    fadir = os.path.join(a.out, "candidates_oriented")
    os.makedirs(fadir, exist_ok=True)
    stats = Counter()

    with open(os.path.join(a.out, "orientation.tsv"), "w") as rep:
        rep.write("bag\tlen\tn_window_hits\tfrac_plus\tstrand_applied\t"
                  "side\tfrac_modal_side\tflipped\n")
        for p in sorted(glob.glob(os.path.join(a.candidates, "*.fasta"))):
            bag = os.path.basename(p)[:-6]
            seq = read_candidate(p)
            wfa = os.path.join(a.windows, bag, bag + ".fasta")
            wjs = os.path.join(a.windows, bag, bag + ".jsonl")
            if not seq or not os.path.exists(wfa) or not os.path.exists(wjs):
                stats["no_windows"] += 1
                continue
            meta = window_meta(wjs)
            with tempfile.TemporaryDirectory() as td:
                q = os.path.join(td, "q.fa")
                with open(q, "w") as fh:
                    fh.write(">cand\n%s\n" % seq)
                db = os.path.join(td, "db")
                subprocess.run([a.makeblastdb, "-in", wfa, "-dbtype", "nucl",
                                "-out", db], capture_output=True)
                r = subprocess.run(
                    [a.blastn, "-task", "blastn", "-word_size", "11",
                     "-dust", "no", "-query", q, "-db", db,
                     "-outfmt", "6 sseqid sstart send pident length",
                     "-max_target_seqs", "5000"],
                    capture_output=True, text=True)
            plus = minus = 0
            sides = Counter()
            for line in r.stdout.splitlines():
                f = line.split("\t")
                if len(f) < 5 or float(f[3]) < a.min_identity:
                    continue
                ss, se = int(f[1]), int(f[2])
                if se > ss:
                    plus += 1
                else:
                    minus += 1
                ao, al = meta.get(f[0], (0, 0))
                lo, hi = (ss, se) if se > ss else (se, ss)
                a0, a1 = lo - 1 - ao, hi - ao
                sides["upstream" if a1 <= 0 else
                      ("downstream" if a0 >= al else "overlaps_tnp")] += 1
            tot = plus + minus
            if not tot:
                stats["no_hit"] += 1
                rep.write("%s\t%d\t0\t\t\t\t\t\n" % (bag, len(seq)))
                continue
            fp = plus / tot
            flip = fp < 0.5
            out = rc(seq) if flip else seq
            side, n = sides.most_common(1)[0]
            stats["flipped" if flip else "kept"] += 1
            with open(os.path.join(fadir, bag + ".fasta"), "w") as fh:
                fh.write(">%s|candidate_tnp_frame|len=%d|strand=%s|side=%s\n%s\n"
                         % (bag, len(out), "minus_of_stored" if flip else
                            "as_stored", side, out))
            rep.write("%s\t%d\t%d\t%.3f\t%s\t%s\t%.3f\t%d\n"
                      % (bag, len(seq), tot, fp,
                         "revcomp" if flip else "as_stored", side,
                         n / tot, int(flip)))
    sys.stderr.write("%s\n" % dict(stats))


if __name__ == "__main__":
    main()
