#!/bin/bash
#SBATCH --job-name=pn_hpwin
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --array=0-47%24
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --time=12:00:00
#SBATCH --output=logs/pn_hpwin_%A_%a.out
#SBATCH --error=logs/pn_hpwin_%A_%a.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Windows for the 48 remote-homology pool bags that were not in the stratified
# run. Same cap rule as that run -- each bag's own distinct-protein count, with
# --cap-dedup -- because a decoy drawn from a differently-capped pool would not
# be comparable to the target it is challenging.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
E=/global/scratch/users/kh36969/protein_ncrna/expand
BAG=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$E/hpool/pool_bags_needing_windows.txt")
CAP=$(awk -F'\t' -v b="$BAG" '$1==b{print $3}' "$E/bag_nr100.tsv")
CAP=${CAP:-1000000}
echo "bag: $BAG  cap: $CAP"
$PY -u scripts/extract_windows.py --census "$E/hpool_census" \
    --genomes out/real_genomes.tsv --family "$BAG" \
    --out "$E/hpool_windows/$BAG" --workers "${SLURM_CPUS_PER_TASK}" \
    --proteins --max-per-clade "$CAP" --cap-mode stratified --cap-dedup
