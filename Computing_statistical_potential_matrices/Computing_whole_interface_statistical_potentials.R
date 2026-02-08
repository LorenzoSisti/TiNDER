#########################################################################################
# THIS SCRIPT COMPUTES AMINO ACID CONTACT MATRICES AND WHOLE-INTERFACE STATISTICAL POTENTIALS FROM A DIRECTORY OF PDB STRUCTURES. 
# IT ANALYZES ANTIBODY–ANTIGEN INTERACTIONS BY:
# 1) CALCULATING SIDE-CHAIN CENTROIDS AND CONTACT MATRICES
# 2) IDENTIFYING INTERFACE RESIDUES BETWEEN ANTIBODY AND ANTIGEN
# 3) SUMMARIZING CONTACT FREQUENCIES ACROSS ALL STRUCTURES
# 4) COMPUTING ASYMMETRIC AND SYMMETRIC WHOLE-INTERFACE STATISTICAL POTENTIALS
# 5) VISUALIZING THE RESULTS AS HEATMAPS.
#########################################################################################

### Required libraries
library(bio3d)
library(dplyr)
library(future)
library(furrr)
library(purrr)
library(progressr)
library(pheatmap)
library(patchwork)
library(ggplotify)
library(reshape2)
library(tidyr)

# Define the path to a custom function files
source("/path/to/your/custom/functions/files/functions.R")

### Set up parallelization to speed up computation
plan(multisession, workers = parallel::detectCores() - 1)
handlers(global = TRUE)
handlers("rstudio")

### Define directories and global parameters
pdb_dir <- "/path/to/your/non/redundant/pdb/files/directory/"
results_dir <- "/path/to/the/directory/where/you/have/the/contact/matrix/and/where/you/want/to/save/the/potential/matrices/"
dir.create(results_dir, showWarnings = FALSE)

# Distance cutoff (Å) to define contact between side-chains centroids
DistCutoff <- 8.5  

# Amino acids ordered according to the Kyte-Doolittle hydrophobicity scale
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE") 

set.seed(1234)

all_pdbs <- list.files(pdb_dir, pattern = "*.pdb", recursive = TRUE, full.names = TRUE)

### Main processing function for each PDB file
whole_interface_pmf_from_pdb_set <- function(pdb_path) {
  
  file_name <- basename(pdb_path)
  
  tryCatch({
    
    # Retrieve and read the PDB file
    pdb_aus <- read.pdb(pdb_path)  
    chain_HL <- c("H", "L")
    chain_AG <- "A"
    
    # The following step resolves issues due to non-standard residue numbering in antibodies, e.g. insertions like SER100, LYS100A, THR100B, ASN101...
    # See the custom functions file for further documentation
    renumbered_df <- renumber_ab_chains(pdb_aus, pdb_path, log_file = "errors.log")
    if (!renumbered_df$ok) {
      return(list(ok = FALSE, filename = file_name, path = pdb_path, error = renumbered_df$error))
    }
    df_coord_renumbered <- renumbered_df$df_coord_renumbered
    
    # Compute centroid coordinates for each residue (average of all side-chains atom coordinates)
    centroidi_df <- as.data.frame(
      df_coord_renumbered %>%
        group_by(chain, resno, resid) %>%
        summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
    )
    
    res_names <- paste(centroidi_df$resid, centroidi_df$resno, centroidi_df$chain, sep = "_")
    df_coord_resid_xyz <- centroidi_df[, c("resid", "resno", "x", "y", "z")] 
    rownames(df_coord_resid_xyz) <- res_names
    
    DistMat <- as.matrix(dist(df_coord_resid_xyz[, c("x", "y", "z")])) # Compute the pairwise Euclidean distance matrix between residue centroids
    DistMat_bin <- ifelse(DistMat <= DistCutoff, 1, 0) # Binarize the distance matrix: 1 if distance ≤ cutoff, 0 otherwise
    
    # Define logical conditions to separate antigen and antibody residues
    condA <- !(centroidi_df$chain %in% c("H", "L")) # Antigen residues
    condHL <- centroidi_df$chain %in% c("H", "L") # Antibody residues
    
    # Subset the distance matrix to contain only antigen–antibody distances
    Inter_DistMat <- DistMat[condA, condHL]
    Inter_DistMat_Bin <- ifelse(Inter_DistMat <= DistCutoff, 1, 0)
    
    BS_A <- apply(Inter_DistMat_Bin, 1, sum)
    BS_A <- BS_A[BS_A != 0]
    BS_HL <- apply(Inter_DistMat_Bin, 2, sum)
    BS_HL <- BS_HL[BS_HL != 0]
    
    # Get names of residues involved in the interface
    BS_A_names <- names(BS_A)
    BS_HL_names <- names(BS_HL)
    
    # Initialize contact matrices (asymmetric and symmetric)
    contact_matrix_asym <- matrix(0, nrow=20, ncol=20, dimnames=list(amino_acids, amino_acids))
    contact_matrix_sym  <- matrix(0, nrow=20, ncol=20, dimnames=list(amino_acids, amino_acids))
    
    # Initialize interface residue counters
    residue_counts_interface_ab <- setNames(rep(0, length(amino_acids)), amino_acids)
    residue_counts_interface_ligando <- setNames(rep(0, length(amino_acids)), amino_acids)
    
    # Populate contact matrices:
    
    # By convention, each contact is counted as bidirectional (i.e., counted twice: once for each residue involved)
    
    #  1) contact_matrix_asym → asymmetric (rows = Ab, columns = Ag)
    #     - preserves the orientation of the contact Ab → Ag
    #     - incremented by +2 for each contact
    
    #  2) contact_matrix_sym → symmetric (Ab–Ag and Ag–Ab are equivalent)
    #     - represents the “classical” symmetric statistical potential logic
    #     - diagonal entries (same amino acid type) are incremented by +2
    #     - off-diagonal pairs (aa1 ≠ aa2) are incremented by +2 in both directions
    
    for (res_ab in BS_HL_names) {
      idx_ab <- which(colnames(Inter_DistMat_Bin) == res_ab)
      contacts <- which(Inter_DistMat_Bin[, idx_ab] == 1)
      for (idx_ag in contacts) {
        res_ag <- rownames(Inter_DistMat_Bin)[idx_ag]
        aa_ab <- strsplit(res_ab, "_")[[1]][1]
        aa_ag <- strsplit(res_ag, "_")[[1]][1]
        if (aa_ab %in% amino_acids && aa_ag %in% amino_acids) {
          
          contact_matrix_asym[aa_ab, aa_ag] <- contact_matrix_asym[aa_ab, aa_ag] + 2
          
          if (aa_ab == aa_ag) {
            contact_matrix_sym[aa_ab, aa_ag] <- contact_matrix_sym[aa_ab, aa_ag] + 2 # homotypic contact
          } else {
            contact_matrix_sym[aa_ab, aa_ag] <- contact_matrix_sym[aa_ab, aa_ag] + 2 # heterotypic contact
            contact_matrix_sym[aa_ag, aa_ab] <- contact_matrix_sym[aa_ag, aa_ab] + 2
          }
        }
      }
    }
    
    
    # Count residue occurrences at the interface
    df_centroidi_anticorpo <- centroidi_df[centroidi_df$chain %in% chain_HL, ]
    rownames(df_centroidi_anticorpo) <- paste(df_centroidi_anticorpo$resid, df_centroidi_anticorpo$resno, df_centroidi_anticorpo$chain, sep = "_")
    df_centroidi_anticorpo_bs <- df_centroidi_anticorpo[rownames(df_centroidi_anticorpo) %in% BS_HL_names, ]
    residue_counts_interface_ab <- residue_counts_interface_ab + table(factor(df_centroidi_anticorpo_bs$resid, levels = amino_acids))
    
    df_centroidi_ligando <- centroidi_df[centroidi_df$chain %in% chain_AG, ]
    rownames(df_centroidi_ligando) <- paste(df_centroidi_ligando$resid, df_centroidi_ligando$resno, df_centroidi_ligando$chain, sep = "_")
    df_centroidi_ligando_bs <- df_centroidi_ligando[rownames(df_centroidi_ligando) %in% BS_A_names, ]
    residue_counts_interface_ligando <- residue_counts_interface_ligando + table(factor(df_centroidi_ligando_bs$resid, levels = amino_acids))
    
    return(list(
      ok = TRUE,
      contact_matrix_asym = contact_matrix_asym,
      contact_matrix_sym  = contact_matrix_sym,
      residue_counts_interface_ab = residue_counts_interface_ab,
      residue_counts_interface_ligando = residue_counts_interface_ligando
    ))
    
  }, error = function(e) {
    return(list(ok = FALSE, error = e$message))
  })
}

### Run the function on all PDB files in parallel
with_progress({
  results_list <- future_map(all_pdbs, whole_interface_pmf_from_pdb_set, .options  = furrr_options(seed = TRUE), .progress = TRUE)
})

### Filter valid results and handle failed files
valid_results <- keep(results_list, ~ .x$ok)
failed_files <- map_chr(discard(results_list, ~ .x$ok), "filename")
cat("Totale file falliti:", length(failed_files), "\n")
if (length(failed_files) > 0) {
  writeLines(failed_files, con = file.path(results_dir, "failed_files.txt"))
}

### Sum all contact matrices and residue counts
contact_matrix_asym_sum <- Reduce(`+`, map(valid_results, "contact_matrix_asym"))
contact_matrix_sym_sum  <- Reduce(`+`, map(valid_results, "contact_matrix_sym"))
residue_counts_interface_ab_sum       <- Reduce(`+`, map(valid_results, "residue_counts_interface_ab"))
residue_counts_interface_ligando_sum  <- Reduce(`+`, map(valid_results, "residue_counts_interface_ligando"))

### Save combined results (CSV and RDS)
datasets <- list(
  contact_matrix_asym_sum              = contact_matrix_asym_sum,
  contact_matrix_sym_sum               = contact_matrix_sym_sum,
  residue_counts_interface_ab_sum      = residue_counts_interface_ab_sum,
  residue_counts_interface_ligando_sum = residue_counts_interface_ligando_sum
)

for (name in names(datasets)) {
  write.csv(datasets[[name]], file = file.path(results_dir, paste0(name, ".csv")), row.names = TRUE)
  saveRDS(datasets[[name]], file = file.path(results_dir, paste0(name, ".rds")))
}

### Load saved data for future sessions ###

#contact_matrix_asym_sum <- readRDS(paste(results_dir, "contact_matrix_asym_sum.rds", sep = ""))
#contact_matrix_sym_sum <- readRDS(paste(results_dir, "contact_matrix_sym_sum.rds", sep = ""))
#residue_counts_interface_ab_sum <- readRDS(file.path(results_dir, "residue_counts_interface_ab_sum.rds"))
#residue_counts_interface_ligando_sum <- readRDS(file.path(results_dir, "residue_counts_interface_ligando_sum.rds"))

### Compute frequencies and statistical potentials

# ASYMMETRIC approach: Compute amino acid frequencies for antibody (paratope) and antigen (epitope)
paratope_freq <- residue_counts_interface_ab_sum / sum(residue_counts_interface_ab_sum)
epitope_freq <- residue_counts_interface_ligando_sum / sum(residue_counts_interface_ligando_sum)

# SYMMETRIC approach: Combine counts from both chains and compute normalized frequencies
par_plus_epi <- residue_counts_interface_ab_sum + residue_counts_interface_ligando_sum
residue_freq <- par_plus_epi / sum(par_plus_epi)

# Compute normalized contact frequencies
asym_lower_tri <- contact_matrix_sym_sum[lower.tri(contact_matrix_sym_sum, diag = TRUE)] # All the symmetrical contact information is in a triangular
contact_freq_sym <- contact_matrix_sym_sum / sum(asym_lower_tri) 
contact_freq_asym <- contact_matrix_asym_sum / sum(contact_matrix_asym_sum)

# Compute statistical potentials
V_asym <- -log(contact_freq_asym / outer(paratope_freq, epitope_freq, "*")) * 2.479
V_sym  <- -log(contact_freq_sym / outer(residue_freq, residue_freq, "*")) * 2.479
V_asym[!is.finite(V_asym)] <- 0
V_sym[!is.finite(V_sym)] <- 0

### Prepare and plot heatmaps
# Compute global min/max across both matrices for consistent color scaling
pot_min <- abs(min(min(V_sym), min(V_asym)))
pot_max <- abs(max(max(V_sym), max(V_asym)))
# Assign row/column names for the asymmetric heatmap
rownames(V_asym) <- paste(amino_acids, "Ab", sep = "_")
colnames(V_asym) <- paste(amino_acids, "Ag", sep = "_")

p_asym <- pheatmap(V_asym,
                   color = colorRampPalette(c("gold1", "white", "dodgerblue2"))(50),
                   breaks = seq(- max(pot_min, pot_max), max(pot_min, pot_max), length.out = 51),
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   main = paste("Asymmetric Potential"),
                   display_numbers = FALSE,
                   fontsize = 10)
#silent = TRUE) 
heatmap_asym_plots <- as.ggplot(p_asym$gtable)

p_sym <- pheatmap(V_sym,
                  color = colorRampPalette(c("gold1", "white", "dodgerblue2"))(50),
                  breaks = seq(- max(pot_min, pot_max), max(pot_min, pot_max), length.out = 51),
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  main = paste("Symmetric Potential"),
                  display_numbers = FALSE,
                  fontsize = 10)
#silent = TRUE) 
heatmap_sym_plots <- as.ggplot(p_sym$gtable) 

# Save as .csv
write.csv(V_asym, file = file.path(results_dir, paste0("V_asym.csv")))
write.csv(V_sym,  file = file.path(results_dir, paste0("V_sym.csv")))

# Prepare and export statistical potentials data frames for further analysis

sym_df <- melt(V_sym)
idx <- which(lower.tri(V_sym, diag = TRUE), arr.ind = TRUE)
aa1 <- rownames(V_sym)[idx[,1]]
aa2 <- colnames(V_sym)[idx[,2]]
pair <- paste(pmin(aa1, aa2), pmax(aa1, aa2), sep = "-")
val  <- V_sym[idx]
sym_pairs <- data.frame(pair = pair, value = as.numeric(val), stringsAsFactors = FALSE)
sym_pairs_sep <- sym_pairs %>%
  separate(pair, into = c("aa1", "aa2"), sep = "-")

write.csv(sym_pairs_sep, file = file.path(results_dir, paste0("V_sym_df.csv")))


asym_df <- melt(V_asym) %>%
  rename(
    aa1 = Var1,
    aa2 = Var2
  )

write.csv(asym_df, file = file.path(results_dir, paste0("V_asym_df.csv")))

# Combine both heatmaps into a single plot
(heatmap_asym_plots | heatmap_sym_plots)
