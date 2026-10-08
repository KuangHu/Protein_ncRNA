"""Deterministic protein-diversity sampling.

Two different jobs live here, and they need opposite samples:

* `farthest_point` — pick a subset that *spans* a clade. Use it when the subset
  has to represent the clade's breadth: capping a clade's members before window
  extraction, choosing BLAST references.
* `hash_sample` — pick a subset that is *representative* of a clade. Use it when
  estimating a clade statistic such as median identity. A diverse subset is the
  wrong tool there: it is deliberately enriched for the clade's outliers, so the
  median identity it reports is biased low, and the bias grows with clade size.

What both replace is slicing an input list, which makes the result a function of
MMseqs output order, manifest order or accession order. That is not random, it
is merely arbitrary, and it correlates with exactly the things one wants
spread — species, assembly batch, submission date.
"""

from __future__ import annotations

import hashlib
from collections import defaultdict


def kmers(s: str, k: int = 5) -> set[str]:
    # A sequence shorter than k becomes its own single token rather than an
    # empty set, so two short proteins can still be compared instead of both
    # scoring 0 against everything.
    return {s[i:i + k] for i in range(len(s) - k + 1)} if len(s) >= k else {s}


def jaccard(a: set[str], b: set[str]) -> float:
    if not a or not b:
        return 0.0
    inter = len(a & b)
    return inter / (len(a) + len(b) - inter)


def hash_sample(ids: list[str], cap: int, seed: int = 0) -> list[str]:
    """A deterministic pseudo-random subset, independent of input order.

    Keyed on the id itself, so the same anchor is kept or dropped no matter how
    the list reached this function, and re-running after a corpus update changes
    the sample only where membership changed.
    """
    if cap <= 0 or len(ids) <= cap:
        return sorted(ids)
    key = lambda i: hashlib.sha1(f"{seed}:{i}".encode()).hexdigest()
    return sorted(sorted(ids, key=key)[:cap])


def farthest_point(ids: list[str], seqs: dict[str, str], cap: int,
                   seed: int = 0, prefer: dict[str, tuple] | None = None,
                   k: int = 5) -> list[str]:
    """Up to `cap` ids spread across the set's own 5-mer diversity.

    Greedy farthest-point: repeatedly take the id least similar to everything
    chosen so far. `prefer` supplies an optional sort key (higher is better,
    e.g. window completeness and length) used to choose the first seed and to
    break exact ties; ties beyond that fall back to a seeded hash of the id, so
    the result does not depend on input order.

    Cost is O(len(ids) x cap) Jaccard comparisons. For the clades this is used
    on -- up to ~9,600 anchors capped to 500 -- that is a few million set
    intersections, tens of seconds, and it runs once per clade.
    """
    if cap <= 0 or len(ids) <= cap:
        return sorted(ids)
    ids = sorted(ids)
    tie = {i: hashlib.sha1(f"{seed}:{i}".encode()).hexdigest() for i in ids}
    pref = prefer or {}
    prof = {i: kmers(seqs[i], k) for i in ids if i in seqs}
    # An id with no protein cannot be placed in protein space. Rather than
    # dropping it silently it goes to the back of the queue, after the
    # diversity-selected ones, so the cap is still filled.
    placed = [i for i in ids if i in prof]
    unplaced = [i for i in ids if i not in prof]
    if not placed:
        return hash_sample(ids, cap, seed)

    first = max(placed, key=lambda i: (pref.get(i, ()), tie[i]))
    chosen, taken = [first], {first}
    mind = {i: 1.0 - jaccard(prof[i], prof[first]) for i in placed}
    while len(chosen) < cap and len(chosen) < len(placed):
        rest = [i for i in placed if i not in taken]
        if not rest:
            break
        # Distance first; where several are equally far -- which is the common
        # case once near-identical members are exhausted and every remaining
        # distance is 0 -- fall back to the preference key, then the hash.
        nxt = max(rest, key=lambda i: (mind[i], pref.get(i, ()), tie[i]))
        chosen.append(nxt)
        taken.add(nxt)
        pn = prof[nxt]
        for i in placed:
            d = 1.0 - jaccard(prof[i], pn)
            if d < mind[i]:
                mind[i] = d
    if len(chosen) < cap and unplaced:
        chosen += hash_sample(unplaced, cap - len(chosen), seed)
    return sorted(chosen)


def _assign(ids: list[str], seeds: list[str], prof: dict[str, set[str]],
            tie: dict[str, str]) -> dict[str, list[str]]:
    """Group each id under its most similar seed.

    The tie-break matters more than it looks. In a clade of near-identical
    proteins every Jaccard is equal, and a bare `max()` then returns whichever
    seed comes first in the list -- so every member lands in one group and the
    rest come back empty. The tie-break must therefore be a hash of the *pair*:
    anything built by concatenating a per-id value with the seed still reduces
    to "largest seed wins" and collapses exactly the same way.
    """
    groups: dict[str, list[str]] = {s: [] for s in seeds}
    for i in ids:
        if i not in prof:
            continue
        groups[max(seeds, key=lambda s: (
            jaccard(prof[i], prof.get(s, set())),
            hashlib.sha1(f"{tie[i]}|{s}".encode()).hexdigest()))].append(i)
    return groups


def stratified_sample(ids: list[str], seqs: dict[str, str], cap: int,
                      seed: int = 0, prefer: dict[str, tuple] | None = None,
                      k: int = 5, dedup: bool = False) -> list[str]:
    """Up to `cap` ids spanning the set's diversity, without taking its outliers.

    This is the right sampler for a *cap*, where `farthest_point` is not.
    Keeping the `cap` mutually most-dissimilar members of a 7,877-member clade
    does not sample its breadth, it selects its 7,877th-percentile oddities --
    fragments, misannotations, remote homologs -- and measured on this corpus
    that was worse than the arbitrary slice it replaced.

    Instead: farthest-point seeds to lay out the space, assign every member to
    its nearest seed, then draw round-robin from the groups, most *typical*
    member of each group first. Breadth comes from visiting every group before
    revisiting any; typicality comes from the within-group ordering. The result
    is exactly `cap` ids whenever the set is larger than that.

    `dedup` collapses exact duplicates before the layout. Two identical
    proteins are at distance 0 by definition, so comparing them adds nothing
    while costing the same as any other pair -- and that cost is the whole run
    time, O(cap x n). On the Arm 2 corpus CDS01936 carries 492,568 anchors over
    22,098 distinct proteins, so 95.5% of the comparisons are between copies:
    measured, ~3.7 hours of work for a choice available in ~12 minutes.

    It is a SPEED option and nothing more. Measured against the full layout it
    covers exactly the same number of distinct sequences -- the plain path
    already seeds distinct points first, because duplicates sit at distance 0 --
    but it does NOT return the same ids, with overlap as low as 23 of 50 on a
    synthetic set. So it is off by default: the frozen genomic baseline
    (v1.1-frozen, job 26743719) was produced without it, and turning it on
    silently would make that run unreproducible. Callers that opt in must record
    that they did.

    Either way the duplicates are not discarded: when the distinct sequences
    run out before the budget does, the existing top-up fills the rest from
    whatever is left, so the result is exactly `cap` ids under both paths.
    """
    if cap <= 0 or len(ids) <= cap:
        return sorted(ids)
    ids = sorted(ids)
    tie = {i: hashlib.sha1(f"{seed}:{i}".encode()).hexdigest() for i in ids}
    prof = {i: kmers(seqs[i], k) for i in ids if i in seqs}
    placed = [i for i in ids if i in prof]
    if not placed:
        return hash_sample(ids, cap, seed)
    pref = prefer or {}

    # One representative per distinct sequence, chosen by the same preference
    # that orders everything else, so a truncated window never represents a
    # sequence that a complete one could.
    by_seq: dict[str, list[str]] = {}
    for i in placed:
        by_seq.setdefault(seqs[i], []).append(i)
    copies: dict[str, list[str]] = {}
    reps: list[str] = []
    for s, mem in by_seq.items():
        mem = sorted(mem, reverse=True, key=lambda i: (pref.get(i, ()), tie[i]))
        reps.append(mem[0])
        copies[mem[0]] = mem[1:]
    reps.sort()
    if dedup and len(reps) < len(placed):
        placed = reps

    seeds = farthest_point(placed, seqs, min(cap, len(placed)), seed=seed,
                           prefer=prefer, k=k)
    groups = _assign(placed, seeds, prof, tie)

    ordered = {}
    for s, mem in groups.items():
        ref = hash_sample(mem, 40, seed)        # typicality probe; cap the cost
        ordered[s] = sorted(mem, reverse=True, key=lambda i: (
            pref.get(i, ()),
            sum(jaccard(prof[i], prof[j]) for j in ref),
            tie[i]))

    out, depth = [], 0
    live = [s for s in seeds if ordered.get(s)]
    while len(out) < cap and live:
        for s in live:
            if len(out) >= cap:
                break
            if depth < len(ordered[s]):
                out.append(ordered[s][depth])
        depth += 1
        live = [s for s in live if depth < len(ordered[s])]
    if len(out) < cap:
        rest = [i for i in ids if i not in set(out)]
        out += hash_sample(rest, cap - len(out), seed)
    return sorted(out[:cap])


def nearest_representatives(ids: list[str], seqs: dict[str, str], k_seeds: int,
                            prefer: dict[str, tuple] | None = None,
                            seed: int = 0, k: int = 5) -> list[str]:
    """Farthest-point seeds, then each seed's group representative.

    Using the farthest-point seeds themselves would make references out of the
    set's oddest members, which is the opposite of what a reference should be.
    """
    if k_seeds <= 0 or len(ids) <= k_seeds:
        return sorted(ids)
    seeds = farthest_point(ids, seqs, k_seeds, seed=seed, prefer=prefer, k=k)
    prof = {i: kmers(seqs[i], k) for i in ids if i in seqs}
    tie = {i: hashlib.sha1(f"{seed}:{i}".encode()).hexdigest() for i in ids}
    pref = prefer or {}
    groups = _assign(ids, seeds, prof, tie)
    out = []
    for s, mem in groups.items():
        if not mem:
            continue
        sample = hash_sample(mem, 40, seed)  # tie-break only; cap the cost
        out.append(max(mem, key=lambda i: (
            pref.get(i, ()),
            sum(jaccard(prof[i], prof[j]) for j in sample),
            tie[i])))
    return sorted(out)[:k_seeds]
