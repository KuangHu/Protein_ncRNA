#!/bin/bash
#SBATCH --job-name=pn_hpool
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --time=08:00:00
#SBATCH --output=logs/pn_hpool_%j.out
#SBATCH --error=logs/pn_hpool_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Remote-homology sibling pools for the 18 discovery/unknown bags that cleared
# the recalibrated floor without ever facing the hard null.
#
# The existing rep_vs_rep.m8 cannot serve. Arm 2 families were built at >=40%
# identity over >=60% coverage of BOTH sequences -- deliberately strict, so a
# shared domain could not merge two unrelated bags. That threshold excludes
# precisely the remote homologs wanted here: these targets are singletons under
# it, which is why they have no null. This search is run at high sensitivity
# with a permissive e-value and the thresholds are applied afterwards, so the
# pool can be tiered rather than cut once.
#
# Both tools are run. MMseqs profile-extended search is fast and gives coverage
# on both sequences; phmmer gives an independent, HMM-based call on the same
# pairs. A pool member found by only one of the two is recorded as such.
MM=/global/home/users/kh36969/.conda/envs/opfi/bin/mmseqs
PHMMER=/global/home/users/kh36969/.conda/envs/opfi/bin/phmmer
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
E=/global/scratch/users/kh36969/protein_ncrna/expand
O=$E/hpool
mkdir -p "$O"
export TMPDIR=/tmp/pnhp_${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== target rep FASTA"
$PY -I - "$E/strata/hpool_targets.txt" "$E/bag_reps.faa" "$O/targets.faa" <<'PYEOF'
import sys
want = set(open(sys.argv[1]).read().split())
keep, n = False, 0
with open(sys.argv[2]) as fh, open(sys.argv[3], "w") as out:
    for line in fh:
        if line.startswith(">"):
            keep = line[1:].split("|")[0] in want
            n += keep
        if keep:
            out.write(line)
print(f"{n} representative proteins for {len(want)} target bags")
PYEOF

echo "=== mmseqs easy-search -s 7.5 (sensitive, permissive e-value)"
$MM easy-search "$O/targets.faa" "$E/bag_reps.faa" "$O/hits_mmseqs.m8" \
    "$TMPDIR/mm" --threads "${SLURM_CPUS_PER_TASK}" \
    -s 7.5 -e 1e-3 --max-seqs 5000 --max-accept 5000 \
    --format-output "query,target,pident,alnlen,evalue,bits,qcov,tcov" \
    >/dev/null
echo "  $(wc -l < "$O/hits_mmseqs.m8") hits"

echo "=== phmmer (independent HMM call on the same corpus)"
$PHMMER --cpu "${SLURM_CPUS_PER_TASK}" -E 1e-3 --domE 1e-3 \
    --tblout "$O/hits_phmmer.tbl" -o /dev/null \
    "$O/targets.faa" "$E/bag_reps.faa"
echo "  $(grep -vc '^#' "$O/hits_phmmer.tbl") hits"

echo "=== build pools"
$PY -u scripts/homology_pool.py --targets "$E/strata/hpool_targets.txt" \
    --mmseqs "$O/hits_mmseqs.m8" --phmmer "$O/hits_phmmer.tbl" \
    --pfam "$E/pfam/bag_pfam.tsv" --nr100 "$E/bag_nr100.tsv" \
    --extracted "$E/strat_windows" --out "$O"
echo "=== outputs: $O"
