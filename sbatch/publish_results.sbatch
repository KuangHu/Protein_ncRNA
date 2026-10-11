#!/bin/bash
#SBATCH --job-name=pn_results
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --time=04:00:00
#SBATCH --output=logs/pn_results_%j.out
#SBATCH --error=logs/pn_results_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Collect everything the two external models produced into results/ so the
# numbers in docs/ can be read, and re-derived, without the scratch tree.
#
# The model score tables are copied as they were written. The benchmark
# tables are RE-RUN here rather than transcribed, so what lands in the repo
# is this run's output and not a remembered number.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
E=/global/scratch/users/kh36969/protein_ncrna/expand
R=results
mkdir -p $R/minerva $R/rnamunin $R/derived $R/benchmarks

# Minerva: eight shards per set, one header. bench = the six-family
# benchmark, broad = the other six families, scan = the open 1-2 protein scan.
for set in bench broad scan; do
  for kind in by_bag regions; do
    out=$R/minerva/${set}_${kind}.tsv
    : > "$out"
    for f in $E/minerva/$set/minerva_${kind}_*.tsv; do
      [ -s "$f" ] || continue
      if [ -s "$out" ]; then tail -n +2 "$f" >> "$out"; else cat "$f" >> "$out"; fi
    done
    echo "$out: $(($(wc -l < "$out") - 1)) rows"
  done
done

# RNAMunin: one record per candidate or cluster member, both strands.
gzip -c "$E/rnamunin/nc_munin_scores.tsv" > $R/rnamunin/nc_munin_scores.tsv.gz
cp "$E/rnamunin/test_scores.tsv" $R/rnamunin/selftest_scores.tsv

# The tables the benchmarks read alongside the two models.
cp "$E/nc_insertion/ir_flag.tsv" $R/derived/ir_flag.tsv
cp "$E/nc_insertion/oriented/orientation.tsv" $R/derived/orientation.tsv
cp "$E/nc_insertion/scan_1_2_protein.tsv" $R/derived/scan_1_2_protein.tsv
cp "$E/nc_insertion/scan_shortlist.tsv" $R/derived/scan_shortlist.tsv

B=$R/benchmarks
$PY -u scripts/separation_bench.py        > $B/separation_bench.txt      2>&1
$PY -u scripts/munin_diagnostics.py       > $B/munin_diagnostics.txt     2>&1
$PY -u scripts/broad_is_bench.py          > $B/broad_is_bench.txt        2>&1
# twice: the sweep's own best cell, then the near-tied median cell the
# write-up quoted, because the two are 0.001 apart on the pooled negative
# and the per-family subsets they give are not
$PY -u scripts/rnamunin_broad_bench.py    > $B/rnamunin_broad_bench.txt  2>&1
$PY -u scripts/rnamunin_broad_bench.py --cell best50,med,0 \
                                          > $B/rnamunin_broad_bench_med.txt 2>&1
$PY -u scripts/final_separation_bench.py  > $B/final_separation_bench.txt 2>&1
wc -l $B/*.txt
du -sh $R
