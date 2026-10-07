#!/bin/bash
#SBATCH --job-name=pn_clades
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=08:00:00
#SBATCH --output=logs/pn_clades_%j.out
#SBATCH --error=logs/pn_clades_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA

PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
CEN=${CEN:-/global/scratch/users/kh36969/protein_ncrna/census_26708202/corrected}
OUT=${OUT:-/global/scratch/users/kh36969/protein_ncrna/clades_v1}
NBR=${NBR:-/global/scratch/users/kh36969/protein_ncrna/nbr_26710278/neighborhood.tsv}
# mmseqs writes many small temp files; keep them node-local.
export TMPDIR=/tmp/pnclades_${SLURM_JOB_ID}
mkdir -p "$TMPDIR" "$OUT"
trap 'rm -rf "$TMPDIR"' EXIT

REUSE=""
[ -f "$OUT/cluster_members.tsv" ] && REUSE="--reuse-clusters"

$PY -u scripts/clade_decompose.py \
    --census "$CEN" --out "$OUT" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --identity --neighborhood "$NBR" $REUSE

$PY -u scripts/selected_clade_summary.py --clades "$OUT"
