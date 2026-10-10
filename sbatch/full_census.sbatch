#!/bin/bash
#SBATCH --job-name=pn_fcens
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --time=04:00:00
#SBATCH --output=logs/pn_fcens_%j.out
#SBATCH --error=logs/pn_fcens_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Census for the whole Arm 2 corpus: all 399 bags. Every family >=4 sibling is
# already inside the corpus, so this list is the corpus itself.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
E=/global/scratch/users/kh36969/protein_ncrna/expand
$PY -u scripts/build_census.py --anchors "$E/anchors" \
    --bags "$E/full/census_bags.txt" --out "$E/full_census"
