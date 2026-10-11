#!/bin/bash
#SBATCH --job-name=pn_ncconf
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --time=12:00:00
#SBATCH --array=0-15
#SBATCH --output=logs/pn_ncconf_%A_%a.out
#SBATCH --error=logs/pn_ncconf_%A_%a.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
export PATH=/global/home/users/kh36969/.conda/envs/opfi/bin:$PATH
E=/global/scratch/users/kh36969/protein_ncrna/expand

# Stage 2. Boundary from the insertions, statistics from the windows: count the
# stage-1 candidate in the UNTRIMMED tnp-anchored windows, 7.1x more loci per
# bag, including the fixed insertions the insertion caller cannot see.
$PY -u scripts/confirm_in_windows.py \
    --candidates "$E/nc_insertion/candidates" \
    --windows "$E/strat_windows" "$E/hpool_windows" "$E/pilot_windows" "$E/full_windows" \
    --out "$E/nc_insertion/window_confirm" \
    --min-identity 80 --min-cov 0.5 --threads 8 \
    --shard "$SLURM_ARRAY_TASK_ID" --nshards 16
