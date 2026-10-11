#!/bin/bash
#SBATCH --job-name=pn_munin
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --account=pc_rubinlab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=06:00:00
#SBATCH --output=/global/scratch/users/kh36969/protein_ncrna/logs/munin_%j.out
#SBATCH --error=/global/scratch/users/kh36969/protein_ncrna/logs/munin_%j.err
set -euo pipefail

ROOT=/global/scratch/users/kh36969/protein_ncrna
TOOLS=$ROOT/tools
OUT=$ROOT/expand/rnamunin
mkdir -p "$OUT" "$ROOT/logs"
export PYTHONNOUSERSITE=1
export TF_CPP_MIN_LOG_LEVEL=3
export OMP_NUM_THREADS=16

VENV=$TOOLS/munin_venv
if [ ! -x "$VENV/bin/python" ]; then
  /global/home/users/kh36969/.conda/envs/opfi/bin/python -m venv "$VENV"
  "$VENV/bin/pip" install -q --no-input --upgrade pip
  "$VENV/bin/pip" install -q --no-input "tensorflow-cpu==2.19.0" "h5py" "numpy<2.2"
fi
"$VENV/bin/python" -c "import keras,tensorflow as tf;print('keras',keras.__version__,'tf',tf.__version__)"

# one FASTA of every stage-1 non-coding candidate and its cluster members
CAND=$ROOT/expand/nc_insertion/candidates
IN=$OUT/nc_candidates.fasta
cat "$CAND"/*.fasta > "$IN"
grep -c '^>' "$IN"

# the authors' own single-sequence path, on one record, as the reference value
head -2 "$IN" > "$OUT/one.fasta"
( cd "$TOOLS/RNAMunin" && "$VENV/bin/python" predict13 "$OUT/one.fasta" ) \
  > "$OUT/one.predict13.tsv"

"$VENV/bin/python" -I /global/home/users/kh36969/tools/Protein_ncRNA/scripts/rnamunin_score.py \
  "$IN" --model "$TOOLS/RNAMunin/Convolution13.keras" \
  --out "$OUT/nc_munin_scores.tsv"

wc -l "$OUT/nc_munin_scores.tsv"
