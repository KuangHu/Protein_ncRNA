#!/usr/bin/env python
"""Flag non-coding candidates whose structure is a terminal inverted repeat.

Minerva scores base pairing, and a transposon end pairs beautifully without
being an RNA: the IS4/IS5 terminal inverted repeats fold into a hairpin as
readily as a ribozyme does. This measures the longest exact inverted repeat
in a candidate -- the arm of the best hairpin -- so a high Minerva score
backed by a 15 bp IR arm can be told apart from one backed by a fold.

This is a REPORTED DIAGNOSTIC, not a filter. Nothing here feeds the score.
"""
import argparse
import glob
import os

COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def longest_ir(s):
    """Longest substring whose reverse complement also occurs in s.

    Arms centred on the sequence's own anti-diagonal are skipped: every
    position trivially matches itself in the reverse complement."""
    s = s.upper()
    rc = s.translate(COMP)[::-1]
    n = len(s)
    best = 0
    prev = [0] * (n + 1)
    for i in range(1, n + 1):
        cur = [0] * (n + 1)
        si = s[i - 1]
        for j in range(1, n + 1):
            if si == rc[j - 1]:
                cur[j] = prev[j - 1] + 1
                if cur[j] > best and abs(i + j - (n + 1)) > 2:
                    best = cur[j]
        prev = cur
    return best


def candidates(cdir):
    for p in sorted(glob.glob(os.path.join(cdir, "*.fasta"))):
        bag = os.path.basename(p)[:-6]
        keep, buf = False, []
        with open(p) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    keep = "|candidate|" in line
                elif keep:
                    buf.append(line)
        if buf:
            yield bag, "".join(buf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-len", type=int, default=600,
                    help="skip longer candidates; the scan is quadratic")
    ap.add_argument("--arm", type=int, default=12,
                    help="arm length at which the hairpin explains the score")
    a = ap.parse_args()
    with open(a.out, "w") as out:
        out.write("bag\tlen\tir_arm\tir_frac\tir_flag\n")
        for bag, s in candidates(a.candidates):
            if len(s) > a.max_len:
                out.write("%s\t%d\t\t\t\n" % (bag, len(s)))
                continue
            r = longest_ir(s)
            out.write("%s\t%d\t%d\t%.4f\t%d\n"
                      % (bag, len(s), r, r / len(s), int(r >= a.arm)))


if __name__ == "__main__":
    main()
