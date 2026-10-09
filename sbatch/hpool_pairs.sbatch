#!/bin/bash
#SBATCH --job-name=pn_hppair
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --time=00:20:00
#SBATCH --output=logs/pn_hppair_%j.out
#SBATCH --error=logs/pn_hppair_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Re-emit the pool with the many-to-many pseudo_family_pairs.tsv. The pool on
# disk was built before that fix and carries only the lossy flat map, which
# would have dropped 46 multi-membership pool bags from every family but one.
# Inputs and flags are byte-identical to hpool_search.sbatch so the pool
# statistics already reported do not move -- only the pairs file is new.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
E=/global/scratch/users/kh36969/protein_ncrna/expand
O=$E/hpool
$PY -u scripts/homology_pool.py --targets "$E/strata/hpool_targets.txt" \
    --mmseqs "$O/hits_mmseqs.m8" --phmmer "$O/hits_phmmer.tbl" \
    --pfam "$E/pfam/bag_pfam.tsv" --nr100 "$E/bag_nr100.tsv" \
    --extracted "$E/strat_windows" --out "$O"
test -s "$O/pseudo_family_pairs.tsv"
echo "pairs: $(( $(wc -l < "$O/pseudo_family_pairs.tsv") - 1 ))"
