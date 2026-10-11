#!/bin/bash
#SBATCH --job-name=pn_orient
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --account=pc_rubinlab
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --time=02:00:00
#SBATCH --output=/global/scratch/users/kh36969/protein_ncrna/logs/orient_%j.out
#SBATCH --error=/global/scratch/users/kh36969/protein_ncrna/logs/orient_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
E=/global/scratch/users/kh36969/protein_ncrna/expand
PYTHONNOUSERSITE=1 /global/home/users/kh36969/.conda/envs/opfi/bin/python -u scripts/orient_candidates.py \
  --candidates "$E/nc_insertion/candidates" \
  --windows "$E/elem_windows" \
  --out "$E/nc_insertion/oriented"
wc -l "$E/nc_insertion/oriented/orientation.tsv"
