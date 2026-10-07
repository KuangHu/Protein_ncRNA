#!/bin/bash
#SBATCH --job-name=pn_rnainst
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --time=03:00:00
#SBATCH --output=logs/pn_rnainst_%j.out
#SBATCH --error=logs/pn_rnainst_%j.err
set -euo pipefail
# R-scape and Infernal go in their own env: opfi already carries the whole
# step 1-4 toolchain and a solver conflict there would strand the pipeline.
__conda_root=${__conda_root:-/global/home/users/kh36969/miniconda3}
. "$__conda_root/etc/profile.d/conda.sh"
conda create -y -n rnafold -c conda-forge -c bioconda rscape infernal mafft viennarna
conda activate rnafold
echo "--- versions ---"
for b in R-scape cmbuild cmsearch cmcalibrate RNAalifold mafft; do
  printf '%-14s ' "$b"; command -v "$b" >/dev/null 2>&1 && ("$b" -h 2>&1 | head -2 | tr '\n' ' ') || echo MISSING
  echo
done
