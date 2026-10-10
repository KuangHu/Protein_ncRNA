#!/bin/bash
#SBATCH --job-name=pn_mtest
#SBATCH --partition=es1
#SBATCH --account=pc_rubinlab
#SBATCH --qos=es_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --gres=gpu:A40:1 --time=02:00:00
#SBATCH --output=logs/pn_mtest_%j.out
#SBATCH --error=logs/pn_mtest_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
PY=/global/home/users/kh36969/.conda/envs/minerva/bin/python
export PYTHONNOUSERSITE=1 HF_HOME=/global/scratch/users/kh36969/hf_cache
E=/global/scratch/users/kh36969/protein_ncrna/expand
$PY -I -u scripts/minerva_score.py --windows "$E/strat_windows" \
    --bags "$E/minerva/test_bags.txt" --blocks "$E/full_run/blocks/blocks.tsv" \
    --out "$E/minerva/test" --max-cases 3
