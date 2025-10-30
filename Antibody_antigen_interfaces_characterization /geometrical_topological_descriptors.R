### THIS SCRIPT EXTRACTS TOPOLOGICAL, GEOMETRIC, AND NUMERIC DESCRIPTORS FROM ANTIBODY–ANTIGEN INTERFACES BASED ON PDB STRUCTURES ###
### OUTPUT IS A SINGLE CSV FILE CONTAINING AVERAGE DESCRIPTORS PER COMPLEX ###

### Libraries ###
library(bio3d)
library(dplyr)
library(future)
library(furrr)
library(purrr)
library(progressr)
library(igraph)

source("/Path/to/your/functions/file/functions.R")

### Set parallelization to speed up computations (use all available CPU but one) ###
plan(multisession, workers = parallel::detectCores() - 1)
handlers(global = TRUE)
handlers("rstudio")

### Setting paths, directories and global variables ###
pdb_dir <- "/Path/to/the/PDB/containing/folder/"
files_dir <- "/Path/where/to/save/outputs/"
DistCutoff <- 8.5  # Two residues are defined as "in contact" if the centroids of their side chain are at less than 8.5 Angstroms

# Amino acids ordered according to the Kyte-Doolittle hydrophobicity scale
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE")

all_pdbs <- list.files(pdb_dir, pattern = "*.pdb", recursive = TRUE, full.names = TRUE)

# Execute
with_progress({
  results_list <- future_map(all_pdbs, geometrical_topological_descriptors, .progress = TRUE, .options = furrr_options(seed = TRUE))
})

valid_results <- results_list[map_lgl(results_list, "ok")]
failed_files <- map_chr(discard(results_list, "ok"), "filename")

# Build a unique dataframe
df_NetDes_general <- bind_rows(map(valid_results, "avg_net_des_df"))

# Save as .csv
write.csv(df_NetDes_general, file = file.path(files_dir, "output_descrittori.csv"), row.names = FALSE)
