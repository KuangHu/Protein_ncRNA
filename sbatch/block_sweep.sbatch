#!/bin/bash
#SBATCH --job-name=pn_sweep
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=04:00:00
#SBATCH --output=logs/pn_sweep_%j.out
#SBATCH --error=logs/pn_sweep_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Recurrence sensitivity: the real-vs-decoy gap must survive the choice, not
# depend on it. Run on the full clade set so the decoy null is the real null.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
WIN=${WIN:-/global/scratch/users/kh36969/protein_ncrna/windows_26714415}
BASE=${BASE:-/global/scratch/users/kh36969/protein_ncrna/sweep_${SLURM_JOB_ID}}
export TMPDIR=/tmp/pnsweep_${SLURM_JOB_ID}
mkdir -p "$TMPDIR" "$BASE"
trap 'rm -rf "$TMPDIR"' EXIT
for R in 0.5 0.3 0.2 0.1; do
  echo "=== min-recurrence $R ==="
  $PY -u scripts/discover_blocks.py --windows "$WIN" --decoys "$WIN/decoys" \
      --out "$BASE/r$R" --workers 8 --threads 4 --tmp "$TMPDIR" \
      --min-recurrence "$R" 2>&1 | tail -22
done
