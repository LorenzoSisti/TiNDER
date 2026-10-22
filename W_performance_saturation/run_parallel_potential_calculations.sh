#!/bin/bash
# Run from the folder containing potential.R and jobs.tsv

# Run w/ nohup ./run_parallel_potential_calculations.sh > run_all.log 2>&1 &

JOBS="/Users/lorenzosisti/TiNDER/data/convergence_asymmetric_GR/jobs.tsv" # Running the layer job first
#JOBS="/Users/lorenzosisti/TiNDER/data/convergence_asymmetric_CDR/jobs.tsv"

tail -n +2 "$JOBS" | cut -f4,5 | while read lst out; do
  echo "Running: $out"
  Rscript /Users/lorenzosisti/TiNDER/knowledge_based_statistical_potentials/computing_g_r_statistical_potential.R "$lst" "$out"
  #Rscript /Users/lorenzosisti/TiNDER/knowledge_based_statistical_potentials/computing_whole_cdr_statistical_potential.R "$lst" "$out"
done
