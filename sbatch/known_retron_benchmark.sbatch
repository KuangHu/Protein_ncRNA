#!/bin/bash
#SBATCH --job-name=pn_retronbench
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=08:00:00
#SBATCH --output=logs/pn_retronbench_%j.out
#SBATCH --error=logs/pn_retronbench_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA

# BENCHMARK BRANCH -- runs only after blind predictions are frozen, and its
# output never feeds selected_clades.tsv.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
CEN=${CEN:-/global/scratch/users/kh36969/protein_ncrna/census_26708202/corrected}
BASE=/global/scratch/users/kh36969/protein_ncrna/known_retron_benchmark_${SLURM_JOB_ID}
mkdir -p "$BASE"

$PY -u -I scripts/import_retrondb.py \
    --csv data/known_retrons/raw/retrondb.csv \
    --out data/known_retrons/known_retrons.tsv

$PY -u scripts/extract_windows.py \
    --census "$CEN" --family RT --out "$BASE/rt_windows" \
    --workers "${SLURM_CPUS_PER_TASK}"

$PY -u scripts/annotate_known_retrons.py \
    --known data/known_retrons/known_retrons.tsv \
    --census "$CEN" --windows "$BASE/rt_windows" \
    --out "$BASE" --threads "${SLURM_CPUS_PER_TASK}"
