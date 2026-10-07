#!/bin/bash
#SBATCH --job-name=pn_dbreport
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=4:00:00
#SBATCH --output=logs/pn_dbreport_%j.out
#SBATCH --error=logs/pn_dbreport_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA

PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
# Write to a separate directory: the census array reads corpus/real_genomes.tsv
# while this runs, and rewriting that file underneath it would mean different
# shards saw different manifests. Promote it deliberately once this finishes.
OUT=/global/scratch/users/kh36969/protein_ncrna/corpus_verified
mkdir -p "$OUT"

# --verify-gzip decompresses every candidate, so this is IO-heavy: SLURM only.
$PY -u scripts/genomes_db_report.py \
    --out "$OUT" --verify-gzip --workers "${SLURM_CPUS_PER_TASK}"

echo "verified manifest rows: $(($(wc -l < "$OUT/real_genomes.tsv") - 1))"
echo
echo "size-filtered manifest (in use by the census):"
echo "  $(($(wc -l < /global/scratch/users/kh36969/protein_ncrna/corpus/real_genomes.tsv) - 1)) rows"
echo "To adopt the verified manifest for future runs:"
echo "  cp $OUT/real_genomes.tsv /global/scratch/users/kh36969/protein_ncrna/corpus/real_genomes.tsv"
