#!/bin/bash
#SBATCH --job-name=pn_mpin
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --time=02:00:00
#SBATCH --output=logs/pn_mpin_%j.out
#SBATCH --error=logs/pn_mpin_%j.err
set -euo pipefail
# --ignore-installed pulled torch cu130 and huggingface-hub 2.x from PyPI
# defaults; minerva-dna pins hub <2, and the cluster's A40/A100 nodes are
# running the cu124 stack the rest of this project uses. Both go back.
ENV=/global/home/users/kh36969/.conda/envs/minerva
export PYTHONNOUSERSITE=1
export PIP_CACHE_DIR=/global/scratch/users/kh36969/pip_cache
export HF_HOME=/global/scratch/users/kh36969/hf_cache
"$ENV/bin/pip" install --force-reinstall \
    torch==2.6.0 --index-url https://download.pytorch.org/whl/cu124
"$ENV/bin/pip" install "huggingface_hub>=0.23,<2" "fsspec<=2026.2.0" \
    pyarrow python-dateutil
"$ENV/bin/pip" check || true
"$ENV/bin/python" -I -c "
import torch, numpy, transformers, huggingface_hub
from minerva.modeling_minerva import MinervaForMaskedLM
from minerva.rna_structure import call_structures
from minerva.sequence_utils import build_prodigal_mixed_sequence
print('torch', torch.__version__, '| hub', huggingface_hub.__version__)
print('imports OK')"
