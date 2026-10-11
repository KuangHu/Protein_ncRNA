#!/bin/bash
#SBATCH --job-name=pn_emptyrb
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --time=24:00:00
#SBATCH --output=logs/pn_emptyrb_%j.out
#SBATCH --error=logs/pn_emptyrb_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
export PATH=/global/home/users/kh36969/.conda/envs/opfi/bin:$PATH
E=/global/scratch/users/kh36969/protein_ncrna/expand

# Widen every empty allele of the 56 target bags from 120 bp to +/-5 kb, so an
# in-element block (median 458 bp from the anchor) actually falls inside it.
$PY -u scripts/rebuild_empty_alleles.py \
    --sites /global/scratch/users/kh36969/fna_ins_discovery/bags_v4/insertions_sites.jsonl \
    --bags data/empty_allele_target_bags.txt \
    --genome-map data/genome_map.tsv \
    --out "$E/empty_wide" \
    --flank 5000 --max-genomes 2000 --max-carriers 5 --threads 16

echo "=== test this run's blocks against the widened empty alleles"
$PY -u scripts/empty_allele_test.py \
    --blocks "$E/elem_run/blocks/blocks.tsv" \
    --block-fasta "$E/elem_run/blocks/blocks" \
    --empty "$E/empty_wide/empty_alleles" \
    --index "$E/empty_wide/empty_allele_index.tsv" \
    --out "$E/empty_wide/test"
