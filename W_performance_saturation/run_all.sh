#!/bin/bash
# Usage: bash run_all.sh [n_rep]
set -u
N_REP=${1:-10}

KB=/Users/lorenzosisti/TiNDER/knowledge_based_statistical_potentials
WP=/Users/lorenzosisti/TiNDER/W_performance_saturation

POT_GR=$KB/computing_g_r_statistical_potential.R                                  # <-- ADATTA (script dei potenziali g(r), con ring_*_potential)
POT_CDR=$KB/computing_whole_cdr_statistical_potential.R
SCORE_GR=$KB/assigning_g_r_W_to_docking_poses.R                          # <-- ADATTA (script che produce layer_score_dockq.csv)
SCORE_CDR=$KB/assigning_cdr_W_to_docking_poses.R                        # <-- ADATTA (script che produce cdr_score_dockq.csv)
MAKE_GROUPS=$WP/make_groups.R
ROC=$WP/roc_batch.R

BASE=/Users/lorenzosisti/TiNDER/data/repeats
DOCKQ=/Users/lorenzosisti/TiNDER/data/dockq_af3.csv

# Stop immediately if a script is missing
for f in "$MAKE_GROUPS" "$POT_GR" "$POT_CDR" "$SCORE_GR" "$SCORE_CDR" "$ROC"; do
  [ -f "$f" ] || { echo "Missing script: $f"; exit 1; }
done

for i in $(seq 1 "$N_REP"); do
  REP=$(printf "%02d" "$i")
  ROOT=$BASE/rep_$REP
  GR=$ROOT/convergence_asymmetric_GR
  CDR=$ROOT/convergence_asymmetric_CDR
  echo "=================== REPETITION $REP ==================="

  # 1. groups (skipped if already created)
  if [ ! -f "$GR/jobs.tsv" ]; then
    Rscript "$MAKE_GROUPS" "$ROOT" "$i" || { echo "make_groups failed"; exit 1; }
  fi

  # 2. potentials, g(r) and CDR, for every group
  tail -n +2 "$GR/jobs.tsv" | cut -f4,5 | while read -r lst out; do
    out_cdr="${out/convergence_asymmetric_GR/convergence_asymmetric_CDR}"
    echo "[POT $REP] $out"
    Rscript "$POT_GR"  "$lst" "$out"     < /dev/null
    Rscript "$POT_CDR" "$lst" "$out_cdr" < /dev/null
  done

  # 3. scoring of the 500 poses + DockQ merge
  while read -r f; do
    echo "[SCORE GR  $REP] $(dirname "$f")"
    Rscript "$SCORE_GR" "$(dirname "$f")" "$DOCKQ" < /dev/null
  done < <(find "$GR" -name ring_sym_potential.csv | sort)

  while read -r f; do
    echo "[SCORE CDR $REP] $(dirname "$f")"
    Rscript "$SCORE_CDR" "$(dirname "$f")" "$DOCKQ" < /dev/null
  done < <(find "$CDR" -name whole_int_sym_potential.csv | sort)

  # 4. ROC AUC per group
  Rscript "$ROC" "$GR"  layer_score_dockq.csv
  Rscript "$ROC" "$CDR" cdr_score_dockq.csv
done

echo "Done."
