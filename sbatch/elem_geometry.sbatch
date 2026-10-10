#!/bin/bash
#SBATCH --job-name=pn_egeo
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --time=08:00:00
#SBATCH --output=logs/pn_egeo_%j.out
#SBATCH --error=logs/pn_egeo_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
export PATH=/global/home/users/kh36969/.conda/envs/opfi/bin:$PATH
E=/global/scratch/users/kh36969/protein_ncrna/expand
$PY -u scripts/element_geometry.py --bags "$E/full/subjects.txt" \
    --windows "$E/strat_windows" "$E/hpool_windows" "$E/pilot_windows" "$E/full_windows" \
    --sites /global/scratch/users/kh36969/fna_ins_discovery/bags_v4/insertions_sites.jsonl \
    --out "$E/element"
