#!/bin/bash
#SBATCH --job-name=pn_a2cens
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --time=04:00:00
#SBATCH --output=logs/pn_a2cens_%j.out
#SBATCH --error=logs/pn_a2cens_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Census for the stratified Arm 2 run: the 117 selected bags plus the 125
# sibling bags of every family >=4 they touch, 242 in all. The siblings are not
# run as subjects; they exist so same_family_permutation has something to
# permute against in each lane.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
E=/global/scratch/users/kh36969/protein_ncrna/expand
$PY -u scripts/build_census.py --anchors "$E/anchors" \
    --bags "$E/strata/census_bags.txt" --out "$E/strat_census"
