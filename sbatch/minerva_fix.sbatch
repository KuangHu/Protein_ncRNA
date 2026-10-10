#!/bin/bash
#SBATCH --job-name=pn_mfix
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --time=02:00:00
#SBATCH --output=logs/pn_mfix_%j.out
#SBATCH --error=logs/pn_mfix_%j.err
set -euo pipefail
# ~/.local/lib/python3.11 shadows this env, so a dependency can look satisfied
# while the env itself lacks it. PYTHONNOUSERSITE + --ignore-installed puts the
# packages where the env can see them, and every run below uses -I so user
# site-packages can never leak back in.
ENV=/global/home/users/kh36969/.conda/envs/minerva
export PYTHONNOUSERSITE=1
export PIP_CACHE_DIR=/global/scratch/users/kh36969/pip_cache
export HF_HOME=/global/scratch/users/kh36969/hf_cache
"$ENV/bin/pip" install --ignore-installed numpy regex scipy accelerate
"$ENV/bin/python" -I -c "
import numpy, regex, transformers, torch
print('numpy', numpy.__version__, '| transformers', transformers.__version__, '| torch', torch.__version__)"
"$ENV/bin/python" -I /global/scratch/users/kh36969/protein_ncrna/expand/minerva/mprobe2.py
