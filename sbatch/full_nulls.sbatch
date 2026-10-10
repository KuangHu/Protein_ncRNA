#!/bin/bash
#SBATCH --job-name=pn_fnull
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --time=24:00:00
#SBATCH --output=logs/pn_fnull_%j.out
#SBATCH --error=logs/pn_fnull_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Assemble the real set, the protein space and both nulls for all 399 subjects.
# Settings are the frozen ones: --replicates 3 --seed 1, permutation_null at
# --min-source-clades 4 --max-source-fraction 0.5, family staging at --min-size
# 4. Nothing here is retuned; the run is wider, not looser.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
E=/global/scratch/users/kh36969/protein_ncrna/expand
P=$E/full_run
REAL=$P/real; DEC=$P/decoys; PROT=$P/proteins
mkdir -p "$REAL" "$DEC" "$PROT"
SUBJ=$E/full/subjects.txt
WDIRS=("$E/strat_windows" "$E/hpool_windows" "$E/pilot_windows" "$E/full_windows")

echo "=== real set: $(wc -l < "$SUBJ") subjects"
MISS=0
while read -r BAG; do
  [ -n "$BAG" ] || continue
  for SRC in "${WDIRS[@]}"; do
    if [ -f "$SRC/$BAG/$BAG.fasta" ]; then
      cp "$SRC/$BAG/$BAG.fasta" "$SRC/$BAG/$BAG.jsonl" "$REAL/"; break
    fi
  done
  [ -f "$REAL/$BAG.fasta" ] || { echo "  MISSING windows: $BAG" >&2; MISS=$((MISS+1)); }
done < "$SUBJ"
echo "  $(ls "$REAL"/*.fasta | wc -l) real clades, $MISS missing"

# ONE protein space for every arm. discover_blocks picks references on
# anchor-protein 5-mer diversity; an arm whose ids miss the lookup falls back to
# DNA and is picking references in a different space from the arm it is compared
# against.
echo "=== anchor protein space"
for SRC in "${WDIRS[@]}"; do
  [ -d "$SRC" ] || continue
  for D in "$SRC"/*/; do
    B=$(basename "$D")
    [ -f "$D/proteins/$B.faa" ] && cp -n "$D/proteins/$B.faa" "$PROT/" || true
  done
done
echo "  $(ls "$PROT" | wc -l) anchor-protein FASTAs"

echo "=== null 1/2: within_bag_shuffle"
while read -r BAG; do
  [ -n "$BAG" ] || continue
  for SRC in "${WDIRS[@]}"; do
    if [ -f "$SRC/$BAG/$BAG.fasta" ]; then
      $PY -u scripts/arm2_nulls.py --windows "$SRC/$BAG" \
          --mode within_bag_shuffle --out "$DEC" --replicates 3 --seed 1 \
          --proteins-in "$SRC/$BAG/proteins" --proteins-out "$PROT"
      break
    fi
  done
done < "$SUBJ"

echo "=== null 2/2: same_family_permutation (family_size >= 4)"
$PY -u scripts/stage_family_windows.py \
    --windows "${WDIRS[@]}" \
    --families "$E/bag_families.tsv" --subjects "$SUBJ" \
    --out "$P/family_staged" --min-size 4
$PY -u scripts/permutation_null.py --windows "$P/family_staged" \
    --out "$P/perm_raw" --replicates 3 --seed 1 \
    --min-source-clades 4 --max-source-fraction 0.5
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
print(f"{n} same_family_permutation decoy clades for the subjects")
PYEOF
cp -n "$P"/family_staged/proteins/*.faa "$PROT/" 2>/dev/null || true
echo "=== real $(ls "$REAL"/*.fasta | wc -l)  decoys $(ls "$DEC"/*.fasta | wc -l)  proteins $(ls "$PROT" | wc -l)"
