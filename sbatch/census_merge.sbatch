#!/bin/bash
#SBATCH --job-name=pn_merge
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=logs/pn_merge_%j.out
#SBATCH --error=logs/pn_merge_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA

PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
OUT=${CENSUS_DIR:?set CENSUS_DIR to the census_<arrayjobid> directory}

$PY -u scripts/merge_shards.py --census "$OUT"
$PY -u scripts/census_report.py --census "$OUT" --threads "${SLURM_CPUS_PER_TASK}"

column -t "$OUT/summary_by_anchor_family.tsv"
echo
column -t "$OUT/divergence/divergence_by_family.tsv"
