################################################################################
# THIS SCRIPT ASSIGNS GLOBAL STATISTICAL POTENTIAL SCORES (SYMMETRIC AND ASYMMETRIC) FOR EACH DOCKING MODEL (PDB FILE). 
# IT ASSIGNS STATISTICAL POTENTIAL VALUES (FROM PRE-COMPUTED MATRICES) TO EACH RESIDUE PAIR FOUND IN CONTACT.
# IT SUMMARIZES AND EXPORTS THE TOTAL AND AVERAGE POTENTIAL VALUES PER MODEL.

# IT WAS USED TO ASSIGN STATISTICAL POTENTIAL VALUES FOR STRCTURES DOCKED WITH ALPHAFOLD3, HDOCK, HADDOCK
################################################################################

### Required libraries
library(bio3d)
library(dplyr)
library(tidyr)
library(future)
library(furrr)

# Define the path to a custom function files
source("/path/to/your/custom/functions/files/functions.R")

### Set up parallelization to speed up computation
plan(multisession, workers = parallel::detectCores() - 1)

### Define directories and global parameters
pdb_dir <- "/path/to/your/docked/structures/files/directory/"
results_dir <- "/path/to/the/directory/where/you/have/the/contact/matrix/and/where/you/want/to/save/the/potential/matrices/"
dir.create(results_dir, showWarnings = FALSE)

# Distance cutoff (Å) to define contact between side-chains centroids
DistCutoff <- 8.5
Nsteps <- 3

# Amino acids ordered according to the Kyte-Doolittle hydrophobicity scale
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE")

# Retrieve all PDB files from the docking folder
all_docked_pdbs <- list.files(pdb_dir, pattern = "*.pdb", recursive = TRUE, full.names = TRUE)

# Load precomputed statistical potential data frames
asym_potentials_whole_int <- read.csv("/path/to/your/asymmetric_potential_df/V_asym_df.csv")
sym_potentials_whole_int <- read.csv("/path/to/your/symmetric_potential_df/V_sym_df.csv")

# Function that computes symmetric statistical potentials for a given PDB model
assign_sym_pmf_to_docking <- function(pdb_path) {
  
  file_name <- basename(pdb_path)
  
  tryCatch({
    
    # Retrieve and read the PDB file
    pdb_aus <- read.pdb(pdb_path)  
    parts <- strsplit(file_name, "_")[[1]]
    chain_HL <- c(parts[2], parts[3])
    chain_AG <- parts[4]
    
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
    #DistMat_bin <- ifelse(DistMat <= DistCutoff, 1, 0) # Binarize the distance matrix: 1 if distance ≤ cutoff, 0 otherwise
    
    # Define logical conditions to separate antigen and antibody residues
    condA <- !(centroidi_df$chain %in% chain_HL) # Antigen residues
    condHL <- centroidi_df$chain %in% chain_HL # Antibody residues
    
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
    
    # Extract contact residue pairs
    contact_indices <- which(Inter_DistMat_Bin == 1, arr.ind = TRUE)
    contacts <- data.frame(
      res1 = rownames(Inter_DistMat_Bin)[contact_indices[, 1]],
      res2 = colnames(Inter_DistMat_Bin)[contact_indices[, 2]]
    )
    
    # Skip if no contacts found
    if (nrow(contacts) == 0) {
      message("⚠️  Nessun contatto trovato in ",pdb_path)
      return(NULL)
    }
    
    # Count how many times each amino acid pair appears in the interface
    contatti_clean <- contacts %>%
      mutate(
        aa1 = sub("_.*", "", res1),
        aa2 = sub("_.*", "", res2)
      ) %>%
      group_by(aa1, aa2) %>%
      summarise(n_contacts = n(), .groups = "drop")
    
    
    # Merge the contatti_clean df with the symmetric potential df and calculate weighted potential
    result <- sym_potentials_whole_int %>%
      select(aa1, aa2, value) %>%
      left_join(contatti_clean, by = c("aa1", "aa2")) %>%
      mutate(
        n_contacts = replace_na(n_contacts, 0),
        potenziale_totale = value * n_contacts
      )
    
    return(result)
    
  }, error = function(e) {
    message("❌ Error in ", pdb_path, ": ", e$message)
    return(NULL)
  })
}

# Function that computes asymmetric statistical potentials (directional) for a PDB model
assign_asym_pmf_to_docking <- function(pdb_path) {
  
  file_name <- basename(pdb_path)
  
  tryCatch({
    
    pdb_aus <- read.pdb(pdb_path)  
    parts <- strsplit(file_name, "_")[[1]]
    chain_HL <- c(parts[2], parts[3])
    chain_AG <- parts[4]
    
    renumbered_df <- renumber_ab_chains(pdb_aus, pdb_path, log_file = "errors.log")
    if (!renumbered_df$ok) {
      return(list(ok = FALSE, filename = file_name, path = pdb_path, error = renumbered_df$error))
    }
    df_coord_renumbered <- renumbered_df$df_coord_renumbered
    
    centroidi_df <- as.data.frame(
      df_coord_renumbered %>%
        group_by(chain, resno, resid) %>%
        summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
    )
    
    res_names <- paste(centroidi_df$resid, centroidi_df$resno, centroidi_df$chain, sep = "_")
    df_coord_resid_xyz <- centroidi_df[, c("resid", "resno", "x", "y", "z")] 
    rownames(df_coord_resid_xyz) <- res_names
    
    DistMat <- as.matrix(dist(df_coord_resid_xyz[, c("x", "y", "z")])) 
    
    condA <- !(centroidi_df$chain %in% chain_HL)
    condHL <- centroidi_df$chain %in% chain_HL
    
    Inter_DistMat <- DistMat[condA, condHL]
    Inter_DistMat_Bin <- ifelse(Inter_DistMat <= DistCutoff, 1, 0)
    
    BS_A <- apply(Inter_DistMat_Bin, 1, sum)
    BS_A <- BS_A[BS_A != 0]
    BS_HL <- apply(Inter_DistMat_Bin, 2, sum)
    BS_HL <- BS_HL[BS_HL != 0]
    
    BS_A_names <- names(BS_A)
    BS_HL_names <- names(BS_HL)
    
    contact_indices <- which(Inter_DistMat_Bin == 1, arr.ind = TRUE)
    contacts <- data.frame(
      res1 = rownames(Inter_DistMat_Bin)[contact_indices[, 1]],
      res2 = colnames(Inter_DistMat_Bin)[contact_indices[, 2]]
    )
    
    if (nrow(contacts) == 0) {
      message("⚠️  Nessun contatto trovato in ",pdb_path)
      return(NULL)
    }
    
    # Annotate each amino acid and to whom it belongs (Ab or Ag)
    contatti_clean <- contacts %>%
      mutate(
        aa1 = paste0(
          sub("_.*", "", res2), "_",
          ifelse(sub(".*_", "", res2) %in% chain_HL, "Ab", "Ag")
        ),
        aa2 = paste0(
          sub("_.*", "", res1), "_",
          ifelse(sub(".*_", "", res1) %in% chain_HL, "Ab", "Ag")
        )
      ) %>%
      group_by(aa1, aa2) %>%
      summarise(n_contacts = n(), .groups = "drop")
    
    result <- asym_potentials_whole_int %>%
      select(aa1, aa2, value) %>%
      left_join(contatti_clean, by = c("aa1", "aa2")) %>%
      mutate(
        n_contacts = replace_na(n_contacts, 0),
        potenziale_totale = value * n_contacts
      )
    
    return(result)
    
  }, error = function(e) {
    message("❌ Error in ", pdb_path, ": ", e$message)
    return(NULL)
  })
}

### Run the function on all docked files in parallel
summary_results <- future_map_dfr(
  all_docked_pdbs,
  function(pdb_path) {
    message("Processing: ", basename(pdb_path))
    
    res_sym  <- assign_sym_pmf_to_docking(pdb_path)
    res_asym <- assign_asym_pmf_to_docking(pdb_path)
    
    # If no contacts were found, return NA values
    if (is.null(res_sym) | is.null(res_asym)) {
      return(data.frame(
        pdb = basename(pdb_path),
        sum_sym = NA,
        sum_asym = NA,
        mean_sym = NA,
        mean_asym = NA
      ))
    }
    
    # Compute total and mean potential values
    data.frame(
      pdb = basename(pdb_path),
      sum_sym  = sum(res_sym$potenziale_totale, na.rm = TRUE),
      sum_asym = sum(res_asym$potenziale_totale, na.rm = TRUE),
      mean_sym  = mean(res_sym$potenziale_totale[res_sym$n_contacts > 0], na.rm = TRUE),
      mean_asym = mean(res_asym$potenziale_totale[res_asym$n_contacts > 0], na.rm = TRUE)
    )
  },
  .progress = TRUE,
  .options  = furrr_options(seed = TRUE)
)

# Save results in .csv
write.csv(summary_results,
          file = file.path(results_dir, "potenziali_per_modello_AF3.csv"),
          row.names = FALSE)

message("\n✅ Analysis completed!")
message("➡️ Results saved to: ", file.path(results_dir, "pmf_for_docked_structures.csv"))
