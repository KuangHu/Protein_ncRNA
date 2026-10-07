#!/bin/bash
#SBATCH --job-name=pn_bench
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=logs/pn_bench_%j.out
#SBATCH --error=logs/pn_bench_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# BENCHMARK BRANCH. Scores frozen blind predictions; writes only under $BENCH.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
BENCH=${BENCH:-/global/scratch/users/kh36969/protein_ncrna/known_retron_benchmark_26714577}
CLADES=${CLADES:-/global/scratch/users/kh36969/protein_ncrna/clades_v1}
BLOCKS=${BLOCKS:-/global/scratch/users/kh36969/protein_ncrna/refsweep_26723281/k3}
WIN=${WIN:-/global/scratch/users/kh36969/protein_ncrna/windows_26714415}
$PY -u scripts/benchmark_clade_overlap.py \
    --bench "$BENCH" --clades "$CLADES" --out "$BENCH" \
    --blocks "$BLOCKS" --windows "$WIN"
