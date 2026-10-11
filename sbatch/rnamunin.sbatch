#!/bin/bash
#SBATCH --job-name=pn_munin
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --account=pc_rubinlab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=06:00:00
#SBATCH --output=/global/scratch/users/kh36969/protein_ncrna/logs/munin_%j.out
#SBATCH --error=/global/scratch/users/kh36969/protein_ncrna/logs/munin_%j.err
set -euo pipefail

ROOT=/global/scratch/users/kh36969/protein_ncrna
TOOLS=$ROOT/tools
OUT=$ROOT/expand/rnamunin
mkdir -p "$OUT" "$ROOT/logs"
export PYTHONNOUSERSITE=1
export PYTHONPATH=$TOOLS/pylibs
export OMP_NUM_THREADS=32
export MKL_NUM_THREADS=32
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python

# every non-coding candidate and cluster member stage 1 emitted
CAND=$ROOT/expand/nc_insertion/candidates
IN=$OUT/nc_candidates.fasta
cat "$CAND"/*.fasta > "$IN"
echo "sequences: $(grep -c '^>' "$IN")"

"$PY" /global/home/users/kh36969/tools/Protein_ncRNA/scripts/rnamunin_score.py \
  "$IN" --model "$TOOLS/RNAMunin/Convolution13.keras" \
  --out "$OUT/nc_munin_scores.tsv" --threads 32 --batch 8192 --selftest

wc -l "$OUT/nc_munin_scores.tsv"
