#!/bin/bash
#SBATCH --job-name=pn_hprun
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --time=48:00:00
#SBATCH --output=logs/pn_hprun_%j.out
#SBATCH --error=logs/pn_hprun_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Hard null for the 18 bags that cleared the recalibrated floor (0.7068) having
# never faced one. Their sibling pool is rebuilt from remote homology rather
# than from the strict Arm 2 family threshold that left them singletons.
#
# Nothing about the scoring changes: identical step-4 settings, identical cap
# rule, and the decoys are produced by the unmodified permutation_null with the
# same --min-source-clades 4 --max-source-fraction 0.5 that Arm 1 used, so its
# strong/weak labels mean the same thing here. They are relabelled
# `homology_pool_permutation` so no table can confuse them with a decoy drawn
# from a real co-recruitment family.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
E=/global/scratch/users/kh36969/protein_ncrna/expand
P=$E/hpool_run
REAL=$P/real; DEC=$P/decoys; BLK=$P/blocks; PROT=$P/proteins
export TMPDIR=/tmp/pnhpr_${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
mkdir -p "$REAL" "$DEC" "$PROT"
T=$E/strata/hpool_targets.txt

echo "=== the 18 targets as the real set"
while read -r BAG; do
  [ -n "$BAG" ] || continue
  for SRC in "$E/strat_windows" "$E/hpool_windows"; do
    if [ -f "$SRC/$BAG/$BAG.fasta" ]; then
      cp "$SRC/$BAG/$BAG.fasta" "$SRC/$BAG/$BAG.jsonl" "$REAL/"; break
    fi
  done
done < "$T"
echo "  $(ls "$REAL"/*.fasta | wc -l) real clades"

echo "=== stage windows as pseudo-families (many-to-many)"
# 46 of 136 pool bags are homologs of more than one target, so a bag appears
# under every pseudo-family it belongs to. --min-size 2 here, not 4: a pool of
# 2 or 3 is a weak null and is labelled weak, but refusing to build it would
# leave the bag with no null at all, which is what this run exists to fix.
$PY -u scripts/stage_family_windows.py \
    --windows "$E/strat_windows" "$E/hpool_windows" \
    --pairs "$E/hpool/pseudo_family_pairs.tsv" --subjects "$T" \
    --out "$P/staged" --min-size 2

echo "=== permutation decoys (Arm 1 criterion, unchanged)"
$PY -u scripts/permutation_null.py --windows "$P/staged" \
    --out "$P/perm_raw" --replicates 3 --seed 1 \
    --min-source-clades 4 --max-source-fraction 0.5

echo "=== keep the targets' decoys, relabel the mode"
$PY -I -u - "$P/perm_raw" "$DEC" "$T" <<'PYEOF'
import shutil, sys
from pathlib import Path
src, dst, tgt = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
want = set(tgt.read_text().split())
n = 0
for p in sorted(src.glob("DECOY_*.fasta")):
    if not any(p.stem.endswith(f"__{b}") for b in want):
        continue
    new = p.stem.replace("DECOY_", "DECOY_homology_pool_permutation__", 1)
    shutil.copy(p, dst / f"{new}.fasta")
    j = p.with_suffix(".jsonl")
    if j.exists():
        shutil.copy(j, dst / f"{new}.jsonl")
    n += 1
print(f"{n} homology_pool_permutation decoy clades for the targets")
PYEOF

echo "=== one protein space for every arm"
for SRC in "$E/strat_windows" "$E/hpool_windows"; do
  for D in "$SRC"/*/; do
    B=$(basename "$D")
    [ -f "$D/proteins/$B.faa" ] && cp -n "$D/proteins/$B.faa" "$PROT/" || true
  done
done
cp -n "$P"/staged/proteins/*.faa "$PROT/" 2>/dev/null || true
echo "  $(ls "$PROT" | wc -l) anchor-protein FASTAs"

echo "=== step 4, identical settings"
$PY -u scripts/discover_blocks.py --windows "$REAL" --decoys "$DEC" \
    --proteins "$PROT" --out "$BLK" --references 3 \
    --min-recurrence 0.5 --workers 8 --threads 4 --tmp "$TMPDIR"

echo "=== null report"
$PY -u scripts/arm2_null_report.py --blocks "$BLK" --decoys "$DEC" \
    --proteins "$PROT" --focus "$T" --out "$P/null_comparison" --floor 0.7068
$PY -u scripts/arm2_floor.py --blocks "$BLK/blocks.tsv" --out "$P/floor"
echo "=== outputs: $P"
