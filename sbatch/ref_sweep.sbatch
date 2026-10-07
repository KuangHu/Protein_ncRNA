#!/bin/bash
#SBATCH --job-name=pn_refsweep
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --time=06:00:00
#SBATCH --output=logs/pn_refsweep_%j.out
#SBATCH --error=logs/pn_refsweep_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# How many references can we add before unrelated decoy members become
# alignable too? min-recurrence fixed at the 0.5 operating point throughout.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
WIN=/global/scratch/users/kh36969/protein_ncrna/windows_26714415
BASE=/global/scratch/users/kh36969/protein_ncrna/refsweep_${SLURM_JOB_ID}
export TMPDIR=/tmp/pnrefs_${SLURM_JOB_ID}
mkdir -p "$TMPDIR" "$BASE"
trap 'rm -rf "$TMPDIR"' EXIT
for K in 3 4 6 8; do
  echo "=== references $K ==="
  $PY -u scripts/discover_blocks.py --windows "$WIN" --decoys "$WIN/decoys" \
      --proteins "$WIN/proteins" --out "$BASE/k$K" --references $K \
      --min-recurrence 0.5 --workers 8 --threads 4 --tmp "$TMPDIR" \
      2>&1 | sed -n '/^{/,/^}/p'
done
