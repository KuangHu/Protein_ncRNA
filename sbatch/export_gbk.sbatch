#!/bin/bash
#SBATCH --job-name=pn_gbk
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --array=0-7%8
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --time=08:00:00
#SBATCH --output=logs/pn_gbk_%A_%a.out
#SBATCH --error=logs/pn_gbk_%A_%a.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# One annotated GenBank per case, one subfolder per bag, for every bag that
# produced a real block. Sharded by bag; each shard streams the census once.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
E=/global/scratch/users/kh36969/protein_ncrna/expand
$PY -I -u scripts/export_loci_gbk.py \
    --windows "$E/strat_windows" "$E/hpool_windows" "$E/pilot_windows" \
    --blocks "$E/strat/blocks/blocks.tsv" \
    --members "$E/strat/blocks/block_members.tsv" \
    --census "$E/strat_census/anchors.tsv" \
    --ranked "$E/hpool_run/ranked/ranked_candidates.tsv" \
    --candidates "$E/strat/report/discovery_candidates.tsv" \
    --pfam "$E/pfam/bag_pfam.tsv" \
    --pfam-sensitive "$E/pfam_sensitive/pfam_sensitive.tsv" \
    --orig /global/scratch/users/kh36969/fna_ins_discovery/bags_v4/insertions_bags.tsv \
    --out "$E/loci_gbk" \
    --shard "${SLURM_ARRAY_TASK_ID}" --nshards "${SLURM_ARRAY_TASK_COUNT}"
