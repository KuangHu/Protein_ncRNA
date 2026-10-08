#!/bin/bash
#SBATCH --job-name=pn_a2win
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --time=08:00:00
#SBATCH --array=0-17
#SBATCH --output=logs/pn_a2win_%A_%a.out
#SBATCH --error=logs/pn_a2win_%A_%a.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Arm 2 null pilot, window extraction. One task per bag: the 6 focus bags plus
# every sibling bag of their families, because the same_family_permutation null
# draws its decoys from those siblings and cannot be built without them.
#
# PILOT CAP. These windows are capped at 2000 anchors, which is NOT the uncapped
# setting chosen for the production run. The pilot's question is which null type
# behaves, and that comparison is only valid if real and decoy sets are the same
# size -- so the cap is held fixed across every arm here and re-opened for
# production once a null is chosen.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
E=/global/scratch/users/kh36969/protein_ncrna/expand
BAG=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$E/pilot_bags.txt")
echo "bag: $BAG"
$PY -u scripts/extract_windows.py --census "$E/pilot_census" \
    --genomes out/real_genomes.tsv --family "$BAG" \
    --out "$E/pilot_windows/$BAG" --workers "${SLURM_CPUS_PER_TASK}" \
    --proteins --max-per-clade "${PILOT_CAP:-2000}" --cap-mode stratified \
    --cap-dedup
