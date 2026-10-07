#!/bin/bash
#SBATCH --job-name=pn_reclass
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --time=01:00:00
#SBATCH --output=logs/pn_reclass_%j.out
#SBATCH --error=logs/pn_reclass_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Relabel only: the verdict column is recomputed from the R-scape numbers
# already in fold_summary.tsv. MAFFT and R-scape are not re-run.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
CAND=${CAND:-/global/scratch/users/kh36969/protein_ncrna/step4_candidates}
WIN=${WIN:-/global/scratch/users/kh36969/protein_ncrna/windows_26714415}
K3=${K3:-/global/scratch/users/kh36969/protein_ncrna/refsweep_26723281/k3}

for F in "$CAND/fold_all080" "$CAND/fold_hc"; do
  [ -f "$F/fold_summary.tsv" ] || continue
  echo "== reclassify $F"
  $PY -u scripts/fold_blocks.py --reclassify "$F/fold_summary.tsv"
done

$PY -u scripts/triage_report.py --candidates "$CAND/step4_candidates.tsv" \
    --blocks "$K3" --windows "$WIN" \
    --fold "$CAND/fold_all080" --padding "$CAND/fasta_all080/block_padding.tsv" \
    --out "$CAND/triage" --merge-window "${MERGE:-500}" --exclude-benchmark-only

# The split IS110 locus cannot be read off a median coordinate; plot every
# member's relative coordinates. Same for the second GroupII_RT clade, which
# shows the flanking architecture independently of the first.
export TMPDIR=/tmp/pnloc_${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
I='IS110__thr0.5__GCA_054628815.1_ENA_JBTPXJ010000055_JBTPXJ010000055.1_54433-55452_-'
G='GroupII_RT__thr0.35__GCA_003903225.1_RUPK01000011.1_48315-50216_-'
$PY -u scripts/locus_report.py --blocks "$K3" --windows "$WIN" \
    --out "$CAND/locus_is110" --fold "$CAND/fold_all080" --tmp "$TMPDIR" \
    --block-id "${I}__B000"
$PY -u scripts/locus_report.py --blocks "$K3" --windows "$WIN" \
    --out "$CAND/locus_groupII_c01" --fold "$CAND/fold_all080" --tmp "$TMPDIR" \
    --block-id "${G}__B000" --block-id "${G}__B001"
