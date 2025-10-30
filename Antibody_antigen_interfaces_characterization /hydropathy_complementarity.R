### SCRIPT TO COMPUTE THE HYDROPATHY COMPLEMENTARITY PHYSICO-CHEMICAL DESCRIPTOR AT THE ANTIBODY-ANTIGEN INTERFACE ###
### OUTPUT INCLUDES A 20×20 CONTACT MATRIX AND AN INTERFACE-WIDE HYDROPATHY SCORE FOR EACH PDB STRUCTURE ###

### Libraries ###
library(bio3d)
library(dplyr)
library(future)
library(furrr)
library(purrr)
library(progressr)
library(ggplot2)

source("/Path/to/your/functions/file/functions.R")

### Set parallelization to speed up computations (use all available CPU but one) ###
plan(multisession, workers = parallel::detectCores() - 1)
handlers(global = TRUE)
handlers("rstudio")

### Setting paths, directories and global variables ###
pdb_dir <- "/path/to/your/non/redundant/pdb/files/directory/"
files_dir <- "/path/to/the/directory/where/you/want/to/save/the/scores/"

DistCutoff <- 8.5  # Two residues are defined as "in contact" if the centroids of their side chain are at less than 8.5 Angstroms

# Amino acids ordered according to the Kyte-Doolittle hydrophobicity scale
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE")

### MAIN ###

Hr_complementarity <- hydropathy_complementarity_matrix(
  milanetti_scale_path = "/Path/to/the/hydrophobicity/scale/MilanettiScale.csv",
  amino_acids = amino_acids
)

Hr_vector    <- Hr_complementarity$Hr_vector
Hr_matrix <- Hr_complementarity$Hr_matrix
print(round(Hr_matrix, 3))

all_pdbs <- list.files(pdb_dir, pattern = "*.pdb", recursive = TRUE, full.names = TRUE)

# Execute
with_progress({
  results_list <- future_map(all_pdbs, avg_hydropathy_complementarity, .progress = TRUE, .options = furrr_options(seed = TRUE))
})

valid_results <- keep(results_list, ~ .x$ok)
failed_files <- map_chr(discard(results_list, ~ .x$ok), "filename")

cat("Total failed files:", length(failed_files), "\n")
if (length(failed_files) > 0) {
  writeLines(failed_files, con = paste0(files_dir, "failed_files.txt"))
}

# Build a data frame with the average interface hydropathy values from successful PDBs
avg_Hr_df <- bind_rows(
  lapply(valid_results, function(res) {
    data.frame(file = res$filename, avg_interface_Hr = res$avg_interface_Hr)
  })
)

# Save as .csv
write.csv(avg_Hr_df, file = paste0(files_dir, "avg_interface_Hr_summary.csv"), row.names = FALSE)

# Plot the distribution of average interface hydropathy values
ggplot(avg_Hr_df, aes(x = avg_interface_Hr)) +
  geom_histogram(binwidth = 0.01, fill = "#0073C2FF", color = "white", alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution of Average Interface Hydropathy (H_r)",
    x = "Average H_r value",
    y = "Number of complexes"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  )
