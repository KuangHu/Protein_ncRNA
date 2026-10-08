#!/bin/bash
#SBATCH --job-name=pn_capdiag
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --time=04:00:00
#SBATCH --output=logs/pn_capdiag_%j.out
#SBATCH --error=logs/pn_capdiag_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
S=/global/scratch/users/kh36969/protein_ncrna
WIN=${WIN:?set WIN to the windows dir}
CEN=${CEN:-$S/census_26708202/corrected}
export TMPDIR=/tmp/pncd_${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
$PY -u scripts/cap_diagnostics.py --windows "$WIN" --census "$CEN" \
    --threads "${SLURM_CPUS_PER_TASK}"
