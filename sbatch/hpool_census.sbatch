#!/bin/bash
#SBATCH --job-name=pn_hpcens
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --time=04:00:00
#SBATCH --output=logs/pn_hpcens_%j.out
#SBATCH --error=logs/pn_hpcens_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
E=/global/scratch/users/kh36969/protein_ncrna/expand
$PY -u scripts/build_census.py --anchors "$E/anchors" \
    --bags "$E/hpool/pool_bags_needing_windows.txt" --out "$E/hpool_census"
