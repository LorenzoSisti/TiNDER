# THIS SCRIPT COMPUTES AMINO-ACID STATISTICAL POTENTIALS WITHIN CONCENTRIC RINGS
# AROUND THE ANTIBODY–ANTIGEN INTERFACE. IT LOADS PRECOMPUTED CONTACT MATRICES
# AND RESIDUE COUNTS, NORMALIZES FREQUENCIES, DERIVES ASYMMETRIC/SYMMETRIC
# POTENTIALS VIA LOG-RATIOS, SAVES CSV OUTPUTS, AND PLOTS HEATMAPS PER RING.

### Libraries ###
library(pheatmap)
library(patchwork)
library(ggplotify)

### Setting paths, directories and global variables ###
pdb_dir <- "/path/to/your/non/redundant/pdb/files/directory/"
results_dir <- "/path/to/the/directory/where/you/want/to/save/the/matrices/"

DistCutoff <- 8.5      # distance cutoff (in Å) to define a contact between residues
Nsteps <- 3            # number of radial rings (shells)

pdbL <- list.files(pdb_dir)

# Amino acids ordered according to the Kyte-Doolittle hydrophobicity scale
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE")

# Initialize storage structures for amino acid frequency vectors per ring
paratope_freq <- vector("list", Nsteps) # Frequencies for paratope residues (antibody)
epitope_freq <- vector("list", Nsteps) # Frequencies for epitope residues (antigen)
par_plus_epi <- vector("list", Nsteps) # Combined counts of both paratope and epitope
residue_freq <- vector("list", Nsteps) # Frequencies for symmetric potential calculation

for (step in 1:Nsteps) {
  paratope_freq[[step]] <- setNames(rep(0, length(amino_acids)), amino_acids)
  epitope_freq[[step]] <- setNames(rep(0, length(amino_acids)), amino_acids)
  par_plus_epi[[step]] <- setNames(rep(0, length(amino_acids)), amino_acids)
  residue_freq[[step]] <- setNames(rep(0, length(amino_acids)), amino_acids)
}

# Initialize table to store the total number of contacts per ring (this should contain the same numbers contained in "count")
emi_sums <- data.frame(
  step = integer(),
  emi_sum = numeric(),
  stringsAsFactors = FALSE
)

# Initialize lists to store heatmap plot objects (ggplot format)
heatmap_asym_plots <- list()
heatmap_sym_plots <- list()

### Load saved data ###
contact_matrix <- readRDS(file.path(results_dir, "name_of_the_contact_matrix.rds"))
residue_counts_interface_ab <- readRDS(file.path(results_dir, "name_of_the_residue_counts_file_ab.rds"))
residue_counts_interface_ligando <- readRDS(file.path(results_dir, "name_of_the_residue_counts_file_ag.rds"))

# Compute statistical potentials for each radial ring
for (step in 1:Nsteps) {
  
  # ASYMMETRIC approach:
  # Compute amino acid frequencies for antibody (paratope) and antigen (epitope)
  paratope_freq[[step]] <- residue_counts_interface_ab[[step]] / sum(residue_counts_interface_ab[[step]])
  epitope_freq[[step]] <- residue_counts_interface_ligando[[step]] / sum(residue_counts_interface_ligando[[step]])
  
  # SYMMETRIC approach:
  # Combine counts from both chains and compute normalized frequencies
  par_plus_epi[[step]] <- residue_counts_interface_ab[[step]] + residue_counts_interface_ligando[[step]]
  residue_freq[[step]] <- par_plus_epi[[step]] / sum(par_plus_epi[[step]])
  
  # Extract and normalize contact matrix for current ring
  contact_matrix_step <- contact_matrix[[step]]
  emi_sum <- sum(contact_matrix_step[lower.tri(contact_matrix_step, diag = TRUE)]) # total number of contacts
  emi_sums <- rbind(emi_sums, data.frame(step = step, emi_sum = emi_sum)) 
  contact_freq <- contact_matrix_step / emi_sum # relative contact frequency
  
  # Compute statistical potentials (asymmetric and symmetric)
  V_asym <- -log(contact_freq / outer(paratope_freq[[step]], epitope_freq[[step]], "*")) * 2.479
  V_sym  <- -log(contact_freq / outer(residue_freq[[step]], residue_freq[[step]], "*")) * 2.479
  
  # Replace infinite/NaN values with 0
  V_asym[!is.finite(V_asym)] <- 0
  V_sym[!is.finite(V_sym)] <- 0
  
  # Compute global min/max across both matrices for consistent color scaling
  pot_min <- abs(min(min(V_sym), min(V_asym)))
  pot_max <- abs(max(max(V_sym), max(V_asym)))
  
  # Assign row/column names for the asymmetric heatmap
  rownames(V_asym) <- paste(amino_acids, "Ab", sep = "_")
  colnames(V_asym) <- paste(amino_acids, "Ag", sep = "_")
  
  # Generate asymmetric heatmap as gtable object (silent mode prevents immediate plotting)
  p_asym <- pheatmap(V_asym,
                     color = colorRampPalette(c("gold1", "white", "dodgerblue2"))(50),
                     breaks = seq(- max(pot_min, pot_max), max(pot_min, pot_max), length.out = 51),
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     main = paste("Asymmetric Potential - Ring", step),
                     display_numbers = FALSE,
                     fontsize = 10,
                     silent = TRUE) 
  heatmap_asym_plots[[step]] <- as.ggplot(p_asym$gtable) # Convert to ggplot object
  
  # Generate symmetric heatmap similarly
  p_sym <- pheatmap(V_sym,
                    color = colorRampPalette(c("gold1", "white", "dodgerblue2"))(50),
                    breaks = seq(- max(pot_min, pot_max), max(pot_min, pot_max), length.out = 51),
                    cluster_rows = FALSE,
                    cluster_cols = FALSE,
                    main = paste("Symmetric Potential - Ring", step),
                    display_numbers = FALSE,
                    fontsize = 10,
                    silent = TRUE) 
  heatmap_sym_plots[[step]] <- as.ggplot(p_sym$gtable) # Convert to ggplot object
  
  # Save as .csv
  write.csv(V_asym, file = file.path(results_dir, paste0("V_asym_step", step, ".csv")))
  write.csv(V_sym,  file = file.path(results_dir, paste0("V_sym_step", step, ".csv")))
}

# Save table of contact sums per ring
write.csv(emi_sums, file = file.path(results_dir, "emi_matrix_sums.csv"), row.names = FALSE)

# Plot all symmetric heatmaps together in a grid layout (3 columns)
combined_sym_plot <- wrap_plots(heatmap_sym_plots, ncol = 3)
print(combined_sym_plot)

# Plot all asymmetric heatmaps together in a grid layout (3 columns)
combined_asym_plot <- wrap_plots(heatmap_asym_plots, ncol = 3) 
print(combined_asym_plot)