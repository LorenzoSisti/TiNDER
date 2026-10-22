#!/bin/bash

# Run w/ nohup ./run_parallel_potential_attribution.sh > run_attribution.log 2>&1 &

SCRIPTS="/Users/lorenzosisti/TiNDER/knowledge_based_statistical_potentials"
GR_ROOT="/Users/lorenzosisti/TiNDER/data/convergence_asymmetric_GR"
CDR_ROOT="/Users/lorenzosisti/TiNDER/data/convergence_asymmetric_CDR"
DOCKQ="/Users/lorenzosisti/TiNDER/data/dockq_af3.csv"

# g(r)
while read -r f; do
  d=$(dirname "$f")
  echo "[GR]  $d"
  Rscript "$SCRIPTS/assigning_g_r_W_to_docking_poses.R" "$d" "$DOCKQ" < /dev/null
done < <(find "$GR_ROOT" -name ring_sym_potential.csv | sort)

# CDR
while read -r f; do
  d=$(dirname "$f")
  echo "[CDR] $d"
  Rscript "$SCRIPTS/assigning_cdr_W_to_docking_poses.R" "$d" "$DOCKQ" < /dev/null
done < <(find "$CDR_ROOT" -name whole_int_sym_potential.csv | sort)
