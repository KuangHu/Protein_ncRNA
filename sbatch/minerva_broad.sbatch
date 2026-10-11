#!/bin/bash
#SBATCH --job-name=pn_mbroad
#SBATCH --partition=es1
#SBATCH --account=pc_rubinlab
#SBATCH --qos=es_normal
#SBATCH --array=0-3%4
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --gres=gpu:A40:1 --time=12:00:00
#SBATCH --output=logs/pn_mbench_%A_%a.out
#SBATCH --error=logs/pn_mbench_%A_%a.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# The six-family benchmark asked three RNA carriers against IS1, IS3 and
# IS4/5/10. This extends the negative to the rest of the IS families present
# in the corpus -- ISL3, IS256, IS66, Mu, Tn3, IS630 -- which already have
# element-trimmed windows but were never folded. Same settings as the bench
# run so the two can be pooled; the labels read the answer, never make it.
#
PY=/global/home/users/kh36969/.conda/envs/minerva/bin/python
export PYTHONNOUSERSITE=1 HF_HOME=/global/scratch/users/kh36969/hf_cache
E=/global/scratch/users/kh36969/protein_ncrna/expand
$PY -I -u scripts/minerva_score.py \
    --windows "$E/elem_windows" \
    --bags data/broad_is_bags.txt \
    --blocks "$E/elem_run/blocks/blocks.tsv" \
    --out "$E/minerva/broad" --max-cases 10 \
    --shard "${SLURM_ARRAY_TASK_ID}" --nshards "${SLURM_ARRAY_TASK_COUNT}"
