#!/bin/bash
#SBATCH --job-name=pn_mver
#SBATCH --partition=es1
#SBATCH --account=pc_rubinlab
#SBATCH --qos=es_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --gres=gpu:A40:1 --time=00:30:00
#SBATCH --output=logs/pn_mver_%j.out
#SBATCH --error=logs/pn_mver_%j.err
set -euo pipefail
PY=/global/home/users/kh36969/.conda/envs/minerva/bin/python
export PYTHONNOUSERSITE=1 HF_HOME=/global/scratch/users/kh36969/hf_cache
$PY -I /global/scratch/users/kh36969/protein_ncrna/expand/minerva/verify.py
