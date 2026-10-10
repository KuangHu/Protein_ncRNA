#!/bin/bash
#SBATCH --job-name=pn_msmoke
#SBATCH --partition=es1
#SBATCH --account=pc_rubinlab
#SBATCH --qos=es_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --gres=gpu:A40:1 --time=02:00:00
#SBATCH --output=logs/pn_msmoke_%j.out
#SBATCH --error=logs/pn_msmoke_%j.err
set -euo pipefail
ENV=/global/home/users/kh36969/.conda/envs/minerva
export PYTHONNOUSERSITE=1 HF_HOME=/global/scratch/users/kh36969/hf_cache
export PIP_CACHE_DIR=/global/scratch/users/kh36969/pip_cache
"$ENV/bin/pip" install -q pandas 2>&1 | tail -2 || true
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
E=/global/scratch/users/kh36969/protein_ncrna/expand
"$ENV/bin/python" -I "$E/minerva/smoke.py" "$E/strat_windows/CDS07978/CDS07978.jsonl"
