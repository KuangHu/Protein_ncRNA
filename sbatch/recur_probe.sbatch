#!/bin/bash
#SBATCH --job-name=pn_rprobe
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --time=24:00:00
#SBATCH --output=logs/pn_rprobe_%j.out
#SBATCH --error=logs/pn_rprobe_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# MEASUREMENT ONLY. This does not change the production pipeline and its output
# is not comparable to any run above -- it exists to put a number on what
# --min-recurrence 0.5 discards, so lowering it can be a decision rather than a
# guess.
#
# 36 scored bags with >=20 distinct anchors produced no block at all at 0.5.
# Either the RNA is not there, or it is there in a minority of the lineages and
# the threshold is hiding it. Those are different results and the run so far
# cannot tell them apart.
#
# Real clades and decoys are the stratified run's own, unchanged, so the only
# difference from that run is the threshold. The decoys come too: a lower
# threshold lets the null produce more blocks as well, and a gain measured
# without that is not a gain.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
E=/global/scratch/users/kh36969/protein_ncrna/expand
P=$E/recur_probe
B=$E/full/recur_probe_bags.txt
export TMPDIR=/tmp/pnrp_${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
rm -rf "$P"; mkdir -p "$P/real" "$P/decoys"

while read -r BAG; do
  [ -n "$BAG" ] || continue
  for EXT in fasta jsonl; do
    [ -f "$E/strat/real/$BAG.$EXT" ] && ln -sf "$E/strat/real/$BAG.$EXT" "$P/real/"
  done
  for F in "$E"/strat/decoys/*__"$BAG".fasta; do
    [ -e "$F" ] || continue
    ln -sf "$F" "$P/decoys/"
    J=${F%.fasta}.jsonl; [ -f "$J" ] && ln -sf "$J" "$P/decoys/"
  done
done < "$B"
echo "probe: $(ls "$P"/real/*.fasta | wc -l) real, $(ls "$P"/decoys/*.fasta 2>/dev/null | wc -l) decoy clades"

for R in 0.25 0.35; do
  echo "=== min-recurrence $R"
  $PY -u scripts/discover_blocks.py --windows "$P/real" --decoys "$P/decoys" \
      --proteins "$E/strat/proteins" --out "$P/r$R" --references 3 \
      --min-recurrence "$R" --workers 8 --threads 4 --tmp "$TMPDIR"
  $PY -u scripts/arm2_floor.py --blocks "$P/r$R/blocks.tsv" --out "$P/r$R/floor"
done
echo "=== outputs: $P"
