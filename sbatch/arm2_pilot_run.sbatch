#!/bin/bash
#SBATCH --job-name=pn_a2pilot
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --time=24:00:00
#SBATCH --output=logs/pn_a2pilot_%j.out
#SBATCH --error=logs/pn_a2pilot_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Arm 2 NULL PILOT. Not a production run and not a result: it asks which of
# three null constructions behaves on this corpus, on six bags chosen to span
# the regimes that matter.
#
# The frozen within-family permutation cannot run on Arm 2's best bags -- 124 of
# 399 are singletons, including the four most diverse, so there is no sibling
# clade to permute against. Three candidate nulls are built and scored side by
# side under identical step-4 settings:
#
#   genomic_background     random same-length windows from the same genomes.
#                          Available for every bag. Easier than sibling
#                          permutation, so it bounds rather than matches.
#   within_bag_shuffle     the bag's own windows with the anchor frame moved
#                          >=2 kb. Strictest: every source of recurrence except
#                          the anchor position survives into the decoy.
#   same_family_permutation  the frozen null, where family_size >= 4 allows it.
#
# The cap is held at the pilot value across real and decoy arms: the comparison
# is only meaningful if the sets are the same size, since recurrence is a
# fraction of members.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
E=/global/scratch/users/kh36969/protein_ncrna/expand
P=$E/pilot
REAL=$P/real
DEC=$P/decoys
BLK=$P/blocks
NT=${SLURM_CPUS_PER_TASK}
REPS=${REPS:-3}
export TMPDIR=/tmp/pna2_${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
mkdir -p "$REAL/proteins" "$DEC"

echo "=== assemble the six focus bags as the real set"
while read -r BAG; do
  [ -n "$BAG" ] || continue
  cp "$E/pilot_windows/$BAG/$BAG.fasta" "$REAL/"
  cp "$E/pilot_windows/$BAG/$BAG.jsonl" "$REAL/"
  cp "$E/pilot_windows/$BAG/proteins/$BAG.faa" "$REAL/proteins/" 2>/dev/null || true
  echo "  $BAG: $(grep -c '^>' "$REAL/$BAG.fasta") windows"
done < "$E/pilot_focus.txt"

echo "=== null 1/3: genomic_background"
while read -r BAG; do
  [ -n "$BAG" ] || continue
  $PY -u scripts/arm2_nulls.py --windows "$E/pilot_windows/$BAG" \
      --mode genomic_background --genomes out/real_genomes.tsv \
      --out "$DEC" --replicates "$REPS" --seed 1
done < "$E/pilot_focus.txt"

echo "=== null 2/3: within_bag_shuffle"
while read -r BAG; do
  [ -n "$BAG" ] || continue
  $PY -u scripts/arm2_nulls.py --windows "$E/pilot_windows/$BAG" \
      --mode within_bag_shuffle --out "$DEC" --replicates "$REPS" --seed 1
done < "$E/pilot_focus.txt"

echo "=== null 3/3: same_family_permutation (family_size >= 4 only)"
$PY -I -u "$E/../_tmp_scripts/stage_family.py"
$PY -u scripts/permutation_null.py --windows "$E/pilot_family_staged" \
    --out "$P/perm_raw" --replicates "$REPS" --seed 1 \
    --min-source-clades 4 --max-source-fraction 0.5
# permutation_null mirrors every staged clade; only the focus bags' decoys are
# wanted here, and they are renamed so the report can tell the three nulls apart
# by clade name alone.
$PY - "$P/perm_raw" "$DEC" "$E/pilot_focus.txt" <<'PYEOF'
import shutil, sys
from pathlib import Path
src, dst, focus = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
want = set(focus.read_text().split())
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
print(f"{n} same_family_permutation decoy clades kept for the focus bags")
PYEOF

echo "=== step 4, identical settings on real + all decoys"
$PY -u scripts/discover_blocks.py --windows "$REAL" --decoys "$DEC" \
    --proteins "$REAL/proteins" --out "$BLK" --references 3 \
    --min-recurrence 0.5 --workers 8 --threads 4 --tmp "$TMPDIR"

echo "=== per-null-type comparison"
$PY -u scripts/arm2_null_report.py --blocks "$BLK" --decoys "$DEC" \
    --focus "$E/pilot_focus.txt" --out "$P/null_comparison"
echo "=== outputs: $P"
