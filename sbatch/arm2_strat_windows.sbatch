#!/bin/bash
#SBATCH --job-name=pn_a2swin
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --array=0-241%24
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --time=12:00:00
#SBATCH --output=logs/pn_a2swin_%A_%a.out
#SBATCH --error=logs/pn_a2swin_%A_%a.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Window extraction for the stratified Arm 2 run, 242 bags.
#
# The cap is each bag's own distinct-protein count, not a flat number. The Arm 2
# corpus is 95.8% duplicate sequence -- 3,710,233 anchors over 154,324 distinct
# proteins across the 117 selected bags -- so a cap set at n_distinct carves one
# window per distinct anchor protein and drops nothing but exact copies. That is
# "no cap on diversity" in the only sense that costs anything: 1.7 Gb of windows
# instead of 40.8 Gb, with 100% of distinct sequences retained in every bag.
#
# What it DOES change: member_fraction and recurrence are then measured over
# distinct proteins rather than raw copies. A block in 1 of 2 distinct anchors
# reads 0.50 here where an uncapped run might read 0.98 because one of the two
# has 40,000 near-identical copies. This run is therefore NOT numerically
# comparable to the v1.1-frozen genomic baseline, and its score floor is
# recalibrated from its own decoys rather than inherited.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
E=/global/scratch/users/kh36969/protein_ncrna/expand

BAG=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$E/strata/census_bags.txt")
# n_distinct for this bag; bags absent from bag_nr100.tsv fall back to a large
# number, which `extract_windows` treats as no cap rather than as zero.
CAP=$(awk -F'\t' -v b="$BAG" '$1==b{print $3}' "$E/bag_nr100.tsv")
CAP=${CAP:-1000000}
echo "bag: $BAG  cap: $CAP (its own distinct-protein count)"

$PY -u scripts/extract_windows.py --census "$E/strat_census" \
    --genomes out/real_genomes.tsv --family "$BAG" \
    --out "$E/strat_windows/$BAG" --workers "${SLURM_CPUS_PER_TASK}" \
    --proteins --max-per-clade "$CAP" --cap-mode stratified --cap-dedup
