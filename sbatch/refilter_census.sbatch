#!/bin/bash
#SBATCH --job-name=pn_refilter
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=logs/pn_refilter_%j.out
#SBATCH --error=logs/pn_refilter_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
CEN=${CEN:?set CEN to the census directory}
$PY -u scripts/refilter_census.py --census "$CEN" --out "$CEN/corrected"
