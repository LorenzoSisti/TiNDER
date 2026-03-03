################################################################################
# THIS SCRIPT ASSIGNS RADIALLY DISTRIBUTED STATISTICAL POTENTIAL SCORES (SYMMETRIC AND ASYMMETRIC) FOR EACH DOCKING MODEL (PDB FILE). 
# IT ASSIGNS STATISTICAL POTENTIAL VALUES (FROM PRE-COMPUTED MATRICES) TO EACH RESIDUE PAIR FOUND IN CONTACT.
# IT SUMMARIZES AND EXPORTS THE TOTAL AND AVERAGE POTENTIAL VALUES PER MODEL.

# IT WAS USED TO ASSIGN STATISTICAL POTENTIAL VALUES FOR STRUCTURES DOCKED WITH ALPHAFOLD3, HDOCK, HADDOCK
################################################################################

### Required libraries
library(bio3d)
library(dplyr)
library(tidyr)
library(future)
library(furrr)

# Define the path to a custom function files
source("/Users/lorenzosisti/Documents/Script_ottimizzati_funzioni/functions.R")

### Set up parallelization to speed up computation
plan(multisession, workers = parallel::detectCores() - 1)

### Define directories and global parameters
#pdb_dir <- "/Users/lorenzosisti/Downloads/AF3_docking/AF3_docking_poses/"
#results_dir <- "/Users/lorenzosisti/Downloads/potenziali_statistici_AF3_prova_strati_opt_marzo_2"
pdb_dir <- "/Users/lorenzosisti/Downloads/models/"
results_dir <- "/Users/lorenzosisti/Downloads/potenziali_statistici_HDOCK_prova_strati_opt_marzo"
dir.create(results_dir, showWarnings = FALSE)

DistCutoff <- 8.5
Nsteps <- 3

# Load precomputed statistical potential data frames
V_asym_wide <- read.csv("/Users/lorenzosisti/Downloads/matrici_stratificate_sparse/V_asym_wide.csv")
V_sym_wide  <- read.csv("/Users/lorenzosisti/Downloads/matrici_stratificate_sparse/V_sym_wide.csv")

# Amino acids ordered according to the Kyte-Doolittle hydrophobicity scale
amino_acids <- c("ARG","LYS","ASN","ASP","GLN","GLU","HIS","PRO","TYR","TRP",
                 "SER","THR","GLY","ALA","MET","CYS","PHE","LEU","VAL","ILE")

# Retrieve all PDB files from the docking folder
all_docked_pdbs <- list.files(pdb_dir, pattern = "*.pdb", recursive = TRUE, full.names = TRUE)

assign_sym_gr_pmf_to_docking <- function(pdb_path) {
  
  file_name <- basename(pdb_path)
  
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
  centroids_df <- as.data.frame(
    df_coord_renumbered %>%
      group_by(chain, resno, resid) %>%
      summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
  )
  
  res_names <- paste(centroids_df$resid, centroids_df$resno, centroids_df$chain, sep = "_")
  df_coord_resid_xyz <- centroids_df[, c("resid", "resno", "x", "y", "z")]
  rownames(df_coord_resid_xyz) <- res_names
  
  DistMat <- as.matrix(dist(df_coord_resid_xyz[, c("x", "y", "z")])) 
  condA <- !(centroids_df$chain %in% chain_HL)
  condHL <- centroids_df$chain %in% chain_HL
  Inter_DistMat <- DistMat[condA, condHL]
  Inter_DistMat_Bin <- ifelse(Inter_DistMat <= DistCutoff, 1, 0)
  
  BS_A <- apply(Inter_DistMat_Bin, 1, sum)
  BS_A <- BS_A[BS_A != 0]
  BS_HL <- apply(Inter_DistMat_Bin, 2, sum)
  BS_HL <- BS_HL[BS_HL != 0]
  if (length(BS_A) == 0 || length(BS_HL) == 0) return(NULL)
  
  BS_A_names <- names(BS_A)
  BS_HL_names <- names(BS_HL)
  df_centroids_BS <- df_coord_resid_xyz[c(BS_A_names, BS_HL_names), ]
  center_BS <- colMeans(df_centroids_BS[, c("x", "y", "z")])
  
  DistMat_centroid <- as.matrix(dist(rbind(df_centroids_BS[, c("x", "y", "z")], center_BS)))
  vet_dist_centroid <- DistMat_centroid[-nrow(DistMat_centroid), ncol(DistMat_centroid)]
  
  borders <- c(0, 0.37, 0.64, 1) * max(vet_dist_centroid)
  contact_indices <- which(Inter_DistMat_Bin == 1, arr.ind = TRUE)
  contacts <- data.frame(
    res1 = rownames(Inter_DistMat_Bin)[contact_indices[, 1]],
    res2 = colnames(Inter_DistMat_Bin)[contact_indices[, 2]]
  )
  contacts$ring_res1 <- findInterval(vet_dist_centroid[contacts$res1], borders, left.open = TRUE)
  contacts$ring_res2 <- findInterval(vet_dist_centroid[contacts$res2], borders, left.open = TRUE)
  
  # --- Calcolo potenziale usando V_sym_wide ---
  contacts <- contacts %>%
    mutate(
      aa1 = sub("_.*", "", res1),
      aa2 = sub("_.*", "", res2),
      pot_col = paste0("V", ring_res1, "_", ring_res2)
    )
  
  # Pesca il valore dalla colonna dinamica di V_sym_wide
  contacts$potenziale_totale <- mapply(function(a1, a2, colname) {
    row_idx <- match(paste(a1, a2, sep="_"), paste(V_sym_wide$aa1, V_sym_wide$aa2, sep="_"))
    if (!is.na(row_idx) && colname %in% colnames(V_sym_wide)) {
      return(V_sym_wide[row_idx, colname])
    } else {
      return(0)
    }
  }, contacts$aa1, contacts$aa2, contacts$pot_col)
  
  total_potential <- sum(contacts$potenziale_totale, na.rm = TRUE)
  n_contacts <- nrow(contacts)
  
  # Calcolo della percentuale di potenziali uguali a zero
  n_zero <- sum(contacts$potenziale_totale == 0, na.rm = TRUE)
  perc_zero <- if (n_contacts > 0) (n_zero / n_contacts) * 100 else NA
  
  #mean_potential <- total_potential / n_contacts
  mean_potential <- if (n_contacts > 0) total_potential / n_contacts else NA
  
  return(list(
    total_potential = total_potential,
    n_contacts = n_contacts,
    mean_potential = mean_potential,
    perc_zero = perc_zero
  ))
}


assign_asym_gr_pmf_to_docking <- function(pdb_path) {
  
  file_name <- basename(pdb_path)
  
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
  centroids_df <- as.data.frame(
    df_coord_renumbered %>%
      group_by(chain, resno, resid) %>%
      summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
  )
  
  res_names <- paste(centroids_df$resid, centroids_df$resno, centroids_df$chain, sep = "_")
  df_coord_resid_xyz <- centroids_df[, c("resid", "resno", "x", "y", "z")]
  rownames(df_coord_resid_xyz) <- res_names
  
  DistMat <- as.matrix(dist(df_coord_resid_xyz[, c("x", "y", "z")])) 
  condA <- !(centroids_df$chain %in% chain_HL)
  condHL <- centroids_df$chain %in% chain_HL
  Inter_DistMat <- DistMat[condA, condHL]
  Inter_DistMat_Bin <- ifelse(Inter_DistMat <= DistCutoff, 1, 0)
  
  BS_A <- apply(Inter_DistMat_Bin, 1, sum)
  BS_A <- BS_A[BS_A != 0]
  BS_HL <- apply(Inter_DistMat_Bin, 2, sum)
  BS_HL <- BS_HL[BS_HL != 0]
  if (length(BS_A) == 0 || length(BS_HL) == 0) return(NULL)
  
  BS_A_names <- names(BS_A)
  BS_HL_names <- names(BS_HL)
  df_centroids_BS <- df_coord_resid_xyz[c(BS_A_names, BS_HL_names), ]
  center_BS <- colMeans(df_centroids_BS[, c("x", "y", "z")])
  
  DistMat_centroid <- as.matrix(dist(rbind(df_centroids_BS[, c("x", "y", "z")], center_BS)))
  vet_dist_centroid <- DistMat_centroid[-nrow(DistMat_centroid), ncol(DistMat_centroid)]
  
  borders <- c(0, 0.37, 0.64, 1) * max(vet_dist_centroid)
  contact_indices <- which(Inter_DistMat_Bin == 1, arr.ind = TRUE)
  contacts <- data.frame(
    res1 = rownames(Inter_DistMat_Bin)[contact_indices[, 1]],
    res2 = colnames(Inter_DistMat_Bin)[contact_indices[, 2]]
  )
  contacts$ring_res1 <- findInterval(vet_dist_centroid[contacts$res1], borders, left.open = TRUE)
  contacts$ring_res2 <- findInterval(vet_dist_centroid[contacts$res2], borders, left.open = TRUE)
  
  # --- Preparazione coppie con suffissi _Ab / _Ag ---
  contacts <- contacts %>%
    mutate(
      aa1 = paste0(sub("_.*", "", res1), "_Ab"),  # chain Ab
      aa2 = paste0(sub("_.*", "", res2), "_Ag"),  # chain Ag
      pot_col = paste0("V", ring_res1, "_", ring_res2)
    )
  
  # --- Pesca potenziale asimmetrico ---
  contacts$potenziale_totale <- mapply(function(a1, a2, colname) {
    row_idx <- match(paste(a1, a2, sep="_"), paste(V_asym_wide$aa1, V_asym_wide$aa2, sep="_"))
    if (!is.na(row_idx) && colname %in% colnames(V_asym_wide)) {
      return(V_asym_wide[row_idx, colname])
    } else {
      return(0)
    }
  }, contacts$aa1, contacts$aa2, contacts$pot_col)
  
  total_potential <- sum(contacts$potenziale_totale, na.rm = TRUE)
  n_contacts <- nrow(contacts)
  
  # Calcolo della percentuale di potenziali uguali a zero
  n_zero <- sum(contacts$potenziale_totale == 0, na.rm = TRUE)
  perc_zero <- if (n_contacts > 0) (n_zero / n_contacts) * 100 else NA
  
  #mean_potential <- total_potential / n_contacts
  mean_potential <- if (n_contacts > 0) total_potential / n_contacts else NA
  
  return(list(
    total_potential = total_potential,
    n_contacts = n_contacts,
    mean_potential = mean_potential,
    perc_zero = perc_zero
  ))
}

summary_results <- future_map_dfr(
  all_docked_pdbs,
  function(pdb_path) {
    message("Processing: ", basename(pdb_path))
    
    res_sym  <- assign_sym_gr_pmf_to_docking(pdb_path)
    res_asym <- assign_asym_gr_pmf_to_docking(pdb_path)
    
    # Se il calcolo fallisce, restituisci NA anche per le nuove colonne
    if (is.null(res_sym) | is.null(res_asym)) {
      return(data.frame(
        pdb = basename(pdb_path),
        sum_sym = NA,
        sum_asym = NA,
        mean_sym = NA,
        mean_asym = NA,
        perc_zero_sym = NA,       # <--- AGGIUNTO
        perc_zero_asym = NA       # <--- AGGIUNTO
      ))
    }
    
    # Se ha successo, estrai il dato calcolato e mettilo nel dataframe
    data.frame(
      pdb = basename(pdb_path),
      sum_sym   = res_sym$total_potential,
      sum_asym  = res_asym$total_potential,
      mean_sym  = res_sym$mean_potential,
      mean_asym = res_asym$mean_potential,
      perc_zero_sym  = res_sym$perc_zero,   # <--- AGGIUNTO
      perc_zero_asym = res_asym$perc_zero   # <--- AGGIUNTO
    )
  },
  .progress = TRUE,
  .options  = furrr_options(seed = TRUE)
)

# Salva in CSV (rimane invariato, il nuovo df avrà in automatico le 2 nuove colonne)
write.csv(summary_results, file = file.path(results_dir, "summary_results_HDOCK_stratificati_sparse.csv"), 
          row.names = FALSE)

# Facoltativo: anche in RDS
saveRDS(summary_results, file = file.path(results_dir, "summary_results_HDOCK_stratificati_sparse.rds"))
