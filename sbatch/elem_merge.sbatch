#!/bin/bash
#SBATCH --job-name=pn_emerge
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --time=08:00:00
#SBATCH --output=logs/pn_emerge_%j.out
#SBATCH --error=logs/pn_emerge_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
export PATH=/global/home/users/kh36969/.conda/envs/opfi/bin:$PATH
E=/global/scratch/users/kh36969/protein_ncrna/expand
P=$E/elem_run
BLK=$P/blocks
NC=$((  $(ls "$P"/real/*.fasta | wc -l) + $(ls "$P"/decoys/*.fasta | wc -l) ))

echo "=== merge $(ls -d "$P"/shard_* | wc -l) shards, expecting $NC clades"
$PY -u scripts/merge_block_shards.py --shards "$P"/shard_*/blocks \
    --out "$BLK" --expect-clades "$NC"

echo "=== recalibrate the floor from this run's own decoys"
# Not inherited. 0.7068 came from 117 bags and 0.7077 from the homology pool;
# a corpus-wide run has its own decoy distribution and assuming either would be
# assuming the answer.
$PY -u scripts/arm2_floor.py --blocks "$BLK/blocks.tsv" --out "$P/floor"
FLOOR=$($PY -I -c "import json;print(json.load(open('$P/floor/score_floor.json'))['floor_no_decoy_support'])")
echo "    floor = $FLOOR"

echo "=== per-null comparison"
$PY -u scripts/arm2_null_report.py --blocks "$BLK" --decoys "$P/decoys" \
    --proteins "$P/proteins" --focus "$E/full/subjects.txt" \
    --out "$P/null_comparison" --floor "$FLOOR"

echo "=== lane tables"
# bag_strata.tsv covers all 399 but carries extra columns; project it to the
# exact shape run_bag_lanes.tsv had so the lane report reads one schema only.
awk -F'\t' 'NR==1{print "bag\tlane\tlane_group\tn_distinct\tfamily_size"; next}
     {print $1"\t"$2"\t"$3"\t"$8"\t"$11}' \
    "$E/strata/bag_strata.tsv" > "$P/bag_lanes.tsv"
echo "    $(( $(wc -l < "$P/bag_lanes.tsv") - 1 )) bags with a lane"
$PY -u scripts/arm2_lane_report.py --blocks "$BLK" --decoys "$P/decoys" \
    --lanes "$P/bag_lanes.tsv" --strata "$E/strata/bag_strata.tsv" \
    --floor "$P/floor/score_floor.json" --out "$P/report"

echo "=== rank every bag that produced a block"
$PY -u scripts/rank_candidates.py \
    --candidates "$P/report/all_bags_annotated.tsv" \
    --null "$P/null_comparison/null_comparison.tsv" \
    --pool "$E/hpool/pool_summary.tsv" --pfam "$E/pfam/bag_pfam.tsv" \
    --pfam-sensitive "$E/pfam_sensitive/pfam_sensitive.tsv" \
    --floor "$FLOOR" --out "$P/ranked"
echo "=== insertion acceptance check on this run's own blocks"
# Expected to accept nearly everything now: the search space IS the element.
# A high host_flank rate here would mean the trim did not work.
$PY -u scripts/insertion_filter.py --blocks "$BLK/blocks.tsv" \
    --block-fasta "$BLK/blocks" \
    --sites /global/scratch/users/kh36969/fna_ins_discovery/bags_v4/insertions_sites.jsonl \
    --pfam "$E/pfam/bag_pfam.tsv" --out "$P/insertion_filter"

echo "=== inverted-repeat check on accepted blocks"
$PY -u scripts/ir_check.py --blocks "$BLK/blocks.tsv" \
    --block-fasta "$BLK/blocks" --pfam "$E/pfam/bag_pfam.tsv" \
    --only-in-insertion "$P/insertion_filter/insertion_filter.tsv" \
    --out "$P/ir_check"

echo "=== outputs: $P"
