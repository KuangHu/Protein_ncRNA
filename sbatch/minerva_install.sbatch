#!/bin/bash
#SBATCH --job-name=pn_minst
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --time=04:00:00
#SBATCH --output=logs/pn_minst_%j.out
#SBATCH --error=logs/pn_minst_%j.err
set -euo pipefail
# A dedicated env. opfi runs the pipeline and must not gain a torch/transformers
# dependency tree mid-analysis; minerva-dna is a step-5 proposer and is kept at
# arm's length from the blind chain in the environment as well as in the code.
export CONDA_PKGS_DIRS=/global/scratch/users/kh36969/conda/pkgs
source /global/home/users/kh36969/.conda/etc/profile.d/conda.sh 2>/dev/null || \
  source "$(conda info --base)/etc/profile.d/conda.sh"
ENV=/global/home/users/kh36969/.conda/envs/minerva
if [ ! -x "$ENV/bin/python" ]; then
  conda create -y -p "$ENV" python=3.11
fi
conda activate "$ENV"
export PIP_CACHE_DIR=/global/scratch/users/kh36969/pip_cache
export HF_HOME=/global/scratch/users/kh36969/hf_cache
pip install --upgrade pip
pip install torch --index-url https://download.pytorch.org/whl/cu124
pip install minerva-dna
python - <<'PY'
import torch, minerva
print("torch", torch.__version__, "cuda", torch.version.cuda)
print("minerva", getattr(minerva, "__version__", "?"))
print("exports:", [x for x in dir(minerva) if not x.startswith("_")][:40])
PY
