#!/bin/bash
#SBATCH --job-name=pn_ir
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --time=04:00:00
#SBATCH --output=logs/pn_ir_%j.out
#SBATCH --error=logs/pn_ir_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
export PATH=/global/home/users/kh36969/.conda/envs/opfi/bin:$PATH
E=/global/scratch/users/kh36969/protein_ncrna/expand
$PY -u scripts/ir_check.py --blocks "$E/full_run/blocks/blocks.tsv" \
    --block-fasta "$E/full_run/blocks/blocks" --pfam "$E/pfam/bag_pfam.tsv" \
    --out "$E/full_run/ir_check" --only-in-insertion "$E/full_run/insertion_filter/insertion_filter.tsv"
