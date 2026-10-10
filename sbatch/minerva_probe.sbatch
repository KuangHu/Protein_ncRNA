#!/bin/bash
#SBATCH --job-name=pn_mprobe
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --time=01:00:00
#SBATCH --output=logs/pn_mprobe_%j.out
#SBATCH --error=logs/pn_mprobe_%j.err
set -euo pipefail
PY=/global/home/users/kh36969/.conda/envs/minerva/bin/python
export HF_HOME=/global/scratch/users/kh36969/hf_cache
$PY -I /global/scratch/users/kh36969/protein_ncrna/expand/minerva/mprobe.py
