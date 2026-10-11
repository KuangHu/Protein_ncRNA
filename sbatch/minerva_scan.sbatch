#!/bin/bash
#SBATCH --job-name=pn_mscan
#SBATCH --partition=es1
#SBATCH --account=pc_rubinlab
#SBATCH --qos=es_normal
#SBATCH --array=0-7%8
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --gres=gpu:A40:1 --time=12:00:00
#SBATCH --output=logs/pn_mbench_%A_%a.out
#SBATCH --error=logs/pn_mbench_%A_%a.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# The open scan. Not a benchmark: every bag whose element encodes one or two
# proteins and that has element-trimmed windows, labelled or not. The known
# RNA families are in here as positive controls -- they are expected near the
# top, and a low rank for them would mean the scan is insensitive, not that
# they are absent. Same settings as bench and broad so all three pool.
#
PY=/global/home/users/kh36969/.conda/envs/minerva/bin/python
export PYTHONNOUSERSITE=1 HF_HOME=/global/scratch/users/kh36969/hf_cache
E=/global/scratch/users/kh36969/protein_ncrna/expand
$PY -I -u scripts/minerva_score.py \
    --windows "$E/elem_windows" \
    --bags data/scan_1_2_protein_bags.txt \
    --blocks "$E/elem_run/blocks/blocks.tsv" \
    --out "$E/minerva/scan" --max-cases 10 \
    --shard "${SLURM_ARRAY_TASK_ID}" --nshards "${SLURM_ARRAY_TASK_COUNT}"
