#!/bin/bash
#SBATCH --job-name=pn_cand
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --time=01:00:00
#SBATCH --output=logs/pn_cand_%j.out
#SBATCH --error=logs/pn_cand_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
BLOCKS=${BLOCKS:-/global/scratch/users/kh36969/protein_ncrna/refsweep_26723281/k3}
BENCH=${BENCH:-/global/scratch/users/kh36969/protein_ncrna/known_retron_benchmark_26714577}
OUT=${OUT:-/global/scratch/users/kh36969/protein_ncrna/step4_candidates}
$PY -u scripts/candidate_table.py --blocks "$BLOCKS" --out "$OUT" --benchmark "$BENCH"
