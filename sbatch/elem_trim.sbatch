#!/bin/bash
#SBATCH --job-name=pn_etrim
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --array=0-7%8
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --time=06:00:00
#SBATCH --output=logs/pn_etrim_%A_%a.out
#SBATCH --error=logs/pn_etrim_%A_%a.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
E=/global/scratch/users/kh36969/protein_ncrna/expand
$PY -u scripts/trim_windows_to_element.py \
    --geometry "$E/element/element_geometry.tsv" \
    --windows "$E/strat_windows" "$E/hpool_windows" "$E/pilot_windows" "$E/full_windows" \
    --out "$E/elem_windows" \
    --shard "${SLURM_ARRAY_TASK_ID}" --nshards "${SLURM_ARRAY_TASK_COUNT}"
