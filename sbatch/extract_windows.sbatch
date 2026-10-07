#!/bin/bash
#SBATCH --job-name=pn_windows
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=06:00:00
#SBATCH --output=logs/pn_windows_%j.out
#SBATCH --error=logs/pn_windows_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA

PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
CEN=${CEN:-/global/scratch/users/kh36969/protein_ncrna/census_26708202/corrected}
CLADES=${CLADES:-/global/scratch/users/kh36969/protein_ncrna/clades_v1}
OUT=${OUT:-/global/scratch/users/kh36969/protein_ncrna/windows_${SLURM_JOB_ID}}
mkdir -p "$OUT"

# No prodigal: anchors.tsv already carries win_start/win_end/strand, so this is
# IO-bound. Cap per clade -- 10k near-identical copies need not all be aligned.
$PY -u scripts/extract_windows.py \
    --census "$CEN" \
    --clades "$CLADES/selected_clades.tsv" \
    --members "$CLADES/cluster_members.tsv" \
    --out "$OUT" --workers "${SLURM_CPUS_PER_TASK}" --proteins \
    --max-per-clade "${MAX_PER_CLADE:-500}"
