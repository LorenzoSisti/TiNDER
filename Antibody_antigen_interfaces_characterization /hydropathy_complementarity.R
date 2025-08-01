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

### Compute the hydrophobicity complementarity 20x20 matrix between each possible pair of residues ###
#### The Milanetti scale was used to retrieve hydropathy values for each amino acid ### 
Hr_scale <- read.csv("/path/to/the/hydrophobicity/scale/MilanettiScale.csv")
Hr_scale$Aa_upper <- toupper(Hr_scale$Aa) # Converti la colonna degli amminoacidi in maiuscolo per il matching
Hr_scale$Aa_upper <- factor(Hr_scale$Aa_upper, levels = amino_acids, ordered = TRUE) # Imposta l'ordine come fattore
Hr_scale_ordered <- Hr_scale[order(Hr_scale$Aa_upper), ] # The Milanetti Scale entries are now ordered according to the Kyte-Doolittle scale
Hr_vector <- setNames(Hr_scale_ordered$Hr, Hr_scale_ordered$Aa_upper) # Create a vector with the Hydrophobicity (Hr) values

# Function to calculate the parabolic function hydropathy complementarity as defined in Grassmann et al. 2025 https://doi.org/10.1021/acs.jcim.4c02286
H_formula <- function(x) {
  -0.033 * (x^2) + 0.363 * x
}

product_matrix <- outer(Hr_vector, Hr_vector, FUN = "*") # Calculate the Hr_vector*Hr_vector product matrix
H_ris_matrix <- H_formula(product_matrix) # Apply the previously defined parabolic function to each matrix element
rownames(H_ris_matrix) <- names(Hr_vector)
colnames(H_ris_matrix) <- names(Hr_vector)
print(round(H_ris_matrix, 3))

pdbL <- list.files(path = pdb_dir, pattern = "*.pdb")

#Function to process each PDB file in the pdbL list and calculate an average hydropathy complementarity value for the whole interface
process_pdb_file <- function(s) {
  tryCatch({
    cat("Elaborating:", pdbL[s], "\n")
    
    # Retrieve and read the PDB file
    path_aus <- paste0(pdb_dir, pdbL[s])
    pdb_aus <- read.pdb(path_aus)
    file_name <- pdbL[s]
    
    # To speed up contact calculations, we use a coarse-grained representation:
    # we retain only the centroids of side chains for each residue.
    # This excludes backbone atoms (N, CA, C, O) for all amino acids,
    # except glycine, where we keep the CA atom (since its side chain is just a hydrogen).
    df_coord <- pdb_aus$atom[pdb_aus$atom$type != "HETATM" & (
      (pdb_aus$atom$resid != "GLY" & !(pdb_aus$atom$elety %in% c("N", "CA", "C", "O"))) |
        (pdb_aus$atom$resid == "GLY" & pdb_aus$atom$elety == "CA")),]
    
    
    if (nrow(df_coord) == 0) {
      writeLines(paste("File whose df_coord is empty:", path_aus), "errors.log", append = TRUE)
      return(list(ok = FALSE, filename = basename(path_aus), error = "df_coord empty"))
    }
    
    # Create a copy of the coordinates to reassign residue numbers per chain
    # This step resolves issues due to non-standard residue numbering in antibodies, e.g. insertions like SER100, LYS100A, THR100B, ASN101...
    df_coord_corrected <- df_coord
    corrected_dfs <- list()
    
    # Reassign residue numbers sequentially within each chain
    for (chain in unique(df_coord_corrected$chain)) {
      residui_unici <- unique(paste(df_coord_corrected$resno[df_coord_corrected$chain == chain], 
                                    df_coord_corrected$insert[df_coord_corrected$chain == chain], sep=""))
      nuova_numerazione <- setNames(seq_along(residui_unici), residui_unici)
      df_coord_corrected$resno[df_coord_corrected$chain == chain] <- nuova_numerazione[paste(df_coord_corrected$resno[df_coord_corrected$chain == chain], 
                                                                                             df_coord_corrected$insert[df_coord_corrected$chain == chain], sep="")]
      corrected_dfs[[chain]] <- df_coord_corrected[df_coord_corrected$chain == chain, ]
    }
    
    df_coord_corrected <- do.call(rbind, corrected_dfs) # Merge corrected data from all chains into a single dataframe
    
    # Compute centroid coordinates for each residue (average of all atom coordinates)
    centroidi_df <- as.data.frame(
      df_coord_corrected %>%
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
    
    # Initialize an amino acid – amino acid contact matrix (20×20)
    contact_matrix <- matrix(0, nrow=20, ncol=20, dimnames=list(amino_acids, amino_acids))
    
    # Loop through each antibody residue in contact
    # By convention, each contact is counted as bidirectional (i.e., counted twice: once for each residue involved)
    for (res_ab in BS_HL_names) {
      idx_ab <- which(colnames(Inter_DistMat_Bin) == res_ab)
      contacts <- which(Inter_DistMat_Bin[, idx_ab] == 1)
      if (length(contacts) > 0) {
        for (idx_ag in contacts) {
          res_ag <- rownames(Inter_DistMat_Bin)[idx_ag]
          aa_ab <- strsplit(res_ab, "_")[[1]][1]
          aa_ag <- strsplit(res_ag, "_")[[1]][1]
          if (aa_ab %in% amino_acids && aa_ag %in% amino_acids) {
            if(aa_ab == aa_ag){
              contact_matrix[aa_ab, aa_ag] <- contact_matrix[aa_ab, aa_ag] + 2 # Increase diagonal entry for homotypic contact
            } else{
              # Update both symmetric entries for heterotypic contact
              contact_matrix[aa_ab, aa_ag] <- contact_matrix[aa_ab, aa_ag] + 2
              contact_matrix[aa_ag, aa_ab] <- contact_matrix[aa_ag, aa_ab] + 2
            }
          }
        }
      }
    }
    
    # Create a lower-triangular version of the contact matrix for further analysis
    emi_matrix <- contact_matrix
    emi_matrix[!lower.tri(contact_matrix, diag = TRUE)] <- 0
    total_contacts <- sum(emi_matrix) # Total number of contacts in the interface
    
    # Apply hydropathy complementarity weights to the contact matrix to get an hydropathy score
    Hydropathy_complementarity_matrix <- contact_matrix * H_ris_matrix
    triangular_matrix <- Hydropathy_complementarity_matrix
    triangular_matrix[!lower.tri(Hydropathy_complementarity_matrix, diag = TRUE)] <- 0
    total_hydropathy <- sum(triangular_matrix)
    avg_interface_Hr <- total_hydropathy/total_contacts
    
    return(list(
      ok = TRUE,
      filename = file_name,
      #filename = pdbL[s],
      avg_interface_Hr = avg_interface_Hr
    ))
    
  }, error = function(e) {
    cat("Error in file:", pdbL[s], "\n", e$message, "\n")
    return(list(ok = FALSE, filename = file_name, error = e$message))
  })
}

# Execute
with_progress({
  results_list <- future_map(1:length(pdbL), process_pdb_file, .progress = TRUE, .options = furrr_options(seed = TRUE))
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
