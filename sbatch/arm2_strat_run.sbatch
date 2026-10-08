#!/bin/bash
#SBATCH --job-name=pn_a2strat
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --time=72:00:00
#SBATCH --output=logs/pn_a2strat_%j.out
#SBATCH --error=logs/pn_a2strat_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Stratified Arm 2 run: 117 subject bags in four lanes, 125 sibling bags staged
# only so the permutation null has a pool.
#
# IS110 / group II intron / IS200-IS605 / TnpB-Cas12 are SPIKE-INS, not the
# target. They are insertion systems already known to carry a structured RNA, so
# what they measure is recovery -- whether this blind caller finds a signal that
# is known to be there. They are reported in their own table and excluded from
# the discovery ranking. The question the run exists to answer is what is left
# once they are removed.
#
# genomic_background is dropped. The pilot scored it in 18 bag-cells and it
# produced zero decoy blocks in all 18: it cannot generate a null distribution,
# and it is the one arm that cannot pick references in protein space. Keeping it
# would only drag any percentile floor toward zero.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
E=/global/scratch/users/kh36969/protein_ncrna/expand
P=$E/strat
REAL=$P/real; DEC=$P/decoys; BLK=$P/blocks; PROT=$P/proteins
REPS=${REPS:-3}
export TMPDIR=/tmp/pna2s_${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
mkdir -p "$REAL" "$DEC" "$PROT"
SUBJ=$E/strata/run_bags.txt

echo "=== assemble the $(wc -l < "$SUBJ") subject bags as the real set"
while read -r BAG; do
  [ -n "$BAG" ] || continue
  if [ -f "$E/strat_windows/$BAG/$BAG.fasta" ]; then
    cp "$E/strat_windows/$BAG/$BAG.fasta" "$E/strat_windows/$BAG/$BAG.jsonl" "$REAL/"
  else
    echo "  MISSING windows for subject $BAG" >&2
  fi
done < "$SUBJ"
echo "  $(ls "$REAL"/*.fasta | wc -l) real clades"

# ONE protein space for every arm, subjects and siblings alike. discover_blocks
# picks references on anchor-protein 5-mer diversity and falls back per id to
# the window's first 3 kb of DNA; an arm whose ids miss the lookup is picking
# references in a different space from the arm it is being compared against.
echo "=== stage anchor proteins for all extracted bags"
for D in "$E"/strat_windows/*/; do
  B=$(basename "$D")
  [ -f "$D/proteins/$B.faa" ] && cp "$D/proteins/$B.faa" "$PROT/" || true
done
echo "  staged $(ls "$PROT" | wc -l) anchor-protein FASTAs"

echo "=== null 1/2: within_bag_shuffle"
while read -r BAG; do
  [ -n "$BAG" ] || continue
  [ -f "$E/strat_windows/$BAG/$BAG.fasta" ] || continue
  $PY -u scripts/arm2_nulls.py --windows "$E/strat_windows/$BAG" \
      --mode within_bag_shuffle --out "$DEC" --replicates "$REPS" --seed 1 \
      --proteins-in "$E/strat_windows/$BAG/proteins" --proteins-out "$PROT"
done < "$SUBJ"

echo "=== null 2/2: same_family_permutation (family_size >= 4 only)"
$PY -u scripts/stage_family_windows.py --windows "$E/strat_windows" \
    --families "$E/bag_families.tsv" --subjects "$SUBJ" \
    --out "$P/family_staged" --min-size 4
$PY -u scripts/permutation_null.py --windows "$P/family_staged" \
    --out "$P/perm_raw" --replicates "$REPS" --seed 1 \
    --min-source-clades 4 --max-source-fraction 0.5
# permutation_null mirrors every staged clade; only the subjects' decoys are
# wanted, renamed so the report can tell the nulls apart by clade name alone.
$PY -I -u - "$P/perm_raw" "$DEC" "$SUBJ" <<'PYEOF'
import shutil, sys
from pathlib import Path
src, dst, subj = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
want = set(subj.read_text().split())
n = 0
for p in sorted(src.glob("DECOY_*.fasta")):
    if not any(p.stem.endswith(f"__{b}") for b in want):
        continue
    new = p.stem.replace("DECOY_", "DECOY_same_family_permutation__", 1)
    shutil.copy(p, dst / f"{new}.fasta")
    j = p.with_suffix(".jsonl")
    if j.exists():
        shutil.copy(j, dst / f"{new}.jsonl")
    n += 1
print(f"{n} same_family_permutation decoy clades kept for the subject bags")
PYEOF

echo "=== step 4, identical settings on real + all decoys"
$PY -u scripts/discover_blocks.py --windows "$REAL" --decoys "$DEC" \
    --proteins "$PROT" --out "$BLK" --references 3 \
    --min-recurrence 0.5 --workers 8 --threads 4 --tmp "$TMPDIR"

echo "=== recalibrate the score floor from this run's own decoys"
# 0.80 came from the genomic corpus at a different cap, where a member fraction
# counted raw copies. Here it counts distinct anchor proteins. Inheriting the
# floor would be assuming the answer.
$PY -u scripts/arm2_floor.py --blocks "$BLK/blocks.tsv" \
    --lanes "$E/strata/run_bag_lanes.tsv" --out "$P/floor"

echo "=== per-null comparison"
$PY -u scripts/arm2_null_report.py --blocks "$BLK" --decoys "$DEC" \
    --proteins "$PROT" --focus "$SUBJ" --out "$P/null_comparison"

echo "=== lane tables"
$PY -u scripts/arm2_lane_report.py --blocks "$BLK" --decoys "$DEC" \
    --lanes "$E/strata/run_bag_lanes.tsv" --strata "$E/strata/bag_strata.tsv" \
    --floor "$P/floor/score_floor.json" --out "$P/report"
echo "=== outputs: $P"
