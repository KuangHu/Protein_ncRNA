#!/bin/bash
#SBATCH --job-name=pn_ncins
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --time=12:00:00
#SBATCH --array=0-15
#SBATCH --output=logs/pn_ncins_%A_%a.out
#SBATCH --error=logs/pn_ncins_%A_%a.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
export PATH=/global/home/users/kh36969/.conda/envs/opfi/bin:$PATH
E=/global/scratch/users/kh36969/protein_ncrna/expand

# insertion = [nc][tnp][nc]. Search the nc regions, nothing else.
$PY -u scripts/nc_rna_in_insertion.py \
    --sites /global/scratch/users/kh36969/fna_ins_discovery/bags_v4/insertions_sites.jsonl \
    --out "$E/nc_insertion" \
    --min-len 50 --min-sites 5 --min-recurrence 0.4 \
    --ident 0.70 --cov 0.70 --decoys 3 --threads 8 \
    --shard "$SLURM_ARRAY_TASK_ID" --nshards 16
