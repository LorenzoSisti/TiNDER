################################################################################
# THIS SCRIPT ASSIGNS RADIALLY DISTRIBUTED STATISTICAL POTENTIAL SCORES (SYMMETRIC AND ASYMMETRIC) FOR EACH DOCKING MODEL (PDB FILE). 
# IT ASSIGNS STATISTICAL POTENTIAL VALUES (FROM PRE-COMPUTED LONG-FORMAT TABLES) TO EACH RESIDUE PAIR FOUND IN CONTACT.
# IT SUMMARIZES AND EXPORTS THE TOTAL AND AVERAGE POTENTIAL VALUES PER MODEL.

# IT WAS USED TO ASSIGN STATISTICAL POTENTIAL VALUES FOR STRUCTURES DOCKED WITH ALPHAFOLD3, HDOCK, HADDOCK
#
# CORREZIONI APPLICATE rispetto alla versione precedente:
# 1) Caricamento dei dati di potenziale aggiornato al formato "long" prodotto dallo
#    script che genera i contatti a ring (ring_asym_potential.csv / ring_sym_potential.csv),
#    al posto dei vecchi file "wide" (V_asym_wide.csv / V_sym_wide.csv).
# 2) assign_sym_gr_pmf_to_docking(): rimossa la logica "_Ab/_Ag" (non necessaria per il
#    potenziale simmetrico, che non distingue la direzione), sostituita con un lookup
#    semplice su aa1/aa2 + ring_pair, usando la tabella simmetrica "specchiata"
#    (df_sym_long_full) per coprire entrambi gli ordinamenti aa1/aa2.
# 3) assign_asym_gr_pmf_to_docking(): corretto il bug di etichettatura per cui res1
#    (che è sempre l'antigene, essendo Inter_DistMat <- DistMat[condA, condHL]) veniva
#    etichettato come "_Ab" e res2 (sempre anticorpo) come "_Ag", cioè al contrario.
#    Ora la chain di appartenenza viene verificata dinamicamente con chain_HL, invece
#    di essere assunta per posizione.
################################################################################

### Required libraries
library(bio3d)
library(dplyr)
library(tidyr)
library(future)
library(furrr)
library(data.table)

# Define the path to a custom function files
source("/Users/lorenzosisti/Documents/Script_ottimizzati_funzioni/functions.R")

### Set up parallelization to speed up computation
plan(multisession, workers = parallel::detectCores() - 1)

### Define directories and global parameters
pdb_dir <- "/Users/lorenzosisti/Downloads/models"
results_dir <- "/Users/lorenzosisti/Downloads/potenziali_statistici_30_06_hdock"
dir.create(results_dir, showWarnings = FALSE)

DistCutoff <- 8.5
Nsteps <- 3

# Amino acids ordered according to the Kyte-Doolittle hydrophobicity scale
amino_acids <- c("ARG","LYS","ASN","ASP","GLN","GLU","HIS","PRO","TYR","TRP",
                 "SER","THR","GLY","ALA","MET","CYS","PHE","LEU","VAL","ILE")

### --- Load precomputed long-format potential tables (from the ring-contacts script) ---
df_asym_long <- fread("/Users/lorenzosisti/Downloads/potenziali_statistici_gr_30_06_data_table_sippl/ring_asym_potential.csv")
df_sym_long  <- fread("/Users/lorenzosisti/Downloads/potenziali_statistici_gr_30_06_data_table_sippl/ring_sym_potential.csv")

# Lookup asimmetrico: chiave diretta resid_ab_resid_ag_ringpair
df_asym_long[, lookup_key := paste(resid_ab, resid_ag, ring_pair, sep = "_")]
setkey(df_asym_long, lookup_key)

# Lookup simmetrico: la tabella ha solo resid_i <= resid_j -> va "specchiata"
# per coprire entrambi gli ordinamenti aa1/aa2 osservati nei contatti
df_sym_long_full <- rbindlist(list(
  df_sym_long[, .(aa_i = resid_i, aa_j = resid_j, ring_pair, potential)],
  df_sym_long[resid_i != resid_j, .(aa_i = resid_j, aa_j = resid_i, ring_pair, potential)]
))
df_sym_long_full[, lookup_key := paste(aa_i, aa_j, ring_pair, sep = "_")]
setkey(df_sym_long_full, lookup_key)

# Retrieve all PDB files from the docking folder
all_docked_pdbs <- list.files(pdb_dir, pattern = "*.pdb", recursive = TRUE, full.names = TRUE)

# pdb_path <- "/Users/lorenzosisti/Downloads/docked_structures_renamed_AF3_11_06/9o58_H_L_A_9.pdb"

### ----------------------------------------------------------------------- ###
### FUNZIONE SIMMETRICA
### ----------------------------------------------------------------------- ###
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
  
  # --- Calcolo potenziale simmetrico (nessuna distinzione di direzione necessaria) ---
  contacts <- contacts %>%
    mutate(
      aa1 = sub("_.*", "", res1), # antigen residue
      aa2 = sub("_.*", "", res2), #antibody residue
      ring_i = pmin(ring_res1, ring_res2),
      ring_j = pmax(ring_res1, ring_res2),
      ring_pair = paste(ring_i, ring_j, sep = "-"),
      lookup_key = paste(aa1, aa2, ring_pair, sep = "_")
    )
  
  contacts$potenziale_totale <- df_sym_long_full[contacts$lookup_key, potential, on = "lookup_key"]
  contacts$potenziale_totale[is.na(contacts$potenziale_totale)] <- 0
  
  total_potential <- sum(contacts$potenziale_totale, na.rm = TRUE)
  n_contacts <- nrow(contacts)
  
  # Calcolo della percentuale di potenziali uguali a zero
  n_zero <- sum(contacts$potenziale_totale == 0, na.rm = TRUE)
  perc_zero <- if (n_contacts > 0) (n_zero / n_contacts) * 100 else NA
  
  mean_potential <- if (n_contacts > 0) total_potential / n_contacts else NA
  
  return(list(
    total_potential = total_potential,
    n_contacts = n_contacts,
    mean_potential = mean_potential,
    perc_zero = perc_zero
  ))
}

### ----------------------------------------------------------------------- ###
### FUNZIONE ASIMMETRICA (corretta)
### ----------------------------------------------------------------------- ###
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
  
  contacts <- contacts %>%
    mutate(
      aa_res1 = sub("_.*", "", res1), # antigen residue
      aa_res2 = sub("_.*", "", res2), # antibody residue
      chain_res1 = sub(".*_", "", res1),
      chain_res2 = sub(".*_", "", res2),
      
      aa_ab = ifelse(chain_res1 %in% chain_HL, aa_res1, aa_res2),  # quello dei due che è in chain_HL -> anticorpo
      aa_ag = ifelse(chain_res1 %in% chain_HL, aa_res2, aa_res1),  # l'altro -> antigene
      
      ring_ab = ifelse(chain_res1 %in% chain_HL, ring_res1, ring_res2),
      ring_ag = ifelse(chain_res1 %in% chain_HL, ring_res2, ring_res1),
      
      ring_i = pmin(ring_ab, ring_ag),
      ring_j = pmax(ring_ab, ring_ag),
      ring_pair = paste(ring_i, ring_j, sep = "-"),
      
      lookup_key = paste(aa_ab, aa_ag, ring_pair, sep = "_")
    )
  
  contacts$potenziale_totale <- df_asym_long[contacts$lookup_key, potential, on = "lookup_key"]
  contacts$potenziale_totale[is.na(contacts$potenziale_totale)] <- 0
  
  total_potential <- sum(contacts$potenziale_totale, na.rm = TRUE)
  n_contacts <- nrow(contacts)
  
  # Calcolo della percentuale di potenziali uguali a zero
  n_zero <- sum(contacts$potenziale_totale == 0, na.rm = TRUE)
  perc_zero <- if (n_contacts > 0) (n_zero / n_contacts) * 100 else NA
  
  mean_potential <- if (n_contacts > 0) total_potential / n_contacts else NA
  
  return(list(
    total_potential = total_potential,
    n_contacts = n_contacts,
    mean_potential = mean_potential,
    perc_zero = perc_zero
  ))
}

### ----------------------------------------------------------------------- ###
### ESECUZIONE IN PARALLELO SU TUTTI I MODELLI DOCKED
### ----------------------------------------------------------------------- ###
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
        perc_zero_sym = NA,
        perc_zero_asym = NA
      ))
    }
    
    # Se ha successo, estrai il dato calcolato e mettilo nel dataframe
    data.frame(
      pdb = basename(pdb_path),
      sum_sym   = res_sym$total_potential,
      sum_asym  = res_asym$total_potential,
      mean_sym  = res_sym$mean_potential,
      mean_asym = res_asym$mean_potential,
      perc_zero_sym  = res_sym$perc_zero,
      perc_zero_asym = res_asym$perc_zero
    )
  },
  .progress = TRUE,
  .options  = furrr_options(seed = TRUE)
)

# Salva in CSV
write.csv(summary_results, file = file.path(results_dir, "potenziali_gr_per_posa.csv"),
          row.names = FALSE)

# Facoltativo: anche in RDS
saveRDS(summary_results, file = file.path(results_dir, "potenziali_gr_per_posa.rds"))
