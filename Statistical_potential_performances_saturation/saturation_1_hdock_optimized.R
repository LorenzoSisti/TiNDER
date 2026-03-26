### BLOCCO 6: SCRIPT DI SCORING PER TUTTE LE POSE (DA LANCIARE CON NOHUP) ###

### --- 1. LIBRERIE E SETUP --- ###
library(bio3d)
library(dplyr)
library(tidyr)
library(future)
library(furrr)
library(purrr)
library(stringr)

# Carica le tue funzioni personalizzate
source("/Users/lorenzosisti/Documents/Script_ottimizzati_funzioni/functions.R")

# Imposta il parallelismo per future_map (usa tutti i core tranne 1)
plan(multisession, workers = parallel::detectCores() - 1)

cat("--- BLOCCO 6: SCRIPT DI SCORING INIZIATO ---\n")
cat("Processori usati:", nbrOfWorkers(), "\n")

### --- 2. DEFINIZIONE DEI PERCORSI (INPUT E OUTPUT) --- ###

# Input 1: La cartella con le LISTE dei potenziali da usare (da Blocco 5)
potentials_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali_marzo/potenziali_statistici_stratificati/"

# Input 2: La cartella con le 22.000 POSE DI DOCKING da analizzare
pdb_dir <- "/Users/lorenzosisti/Downloads/models/"

# Output: La cartella dove salvare i file CSV di risultati
results_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali_marzo/risultati_scoring_finale_HDOCK/"
dir.create(results_dir, showWarnings = FALSE)

# Parametri globali
DistCutoff <- 8.5
Nsteps <- 3
amino_acids <- c("ARG","LYS","ASN","ASP","GLN","GLU","HIS","PRO","TYR","TRP",
                 "SER","THR","GLY","ALA","MET","CYS","PHE","LEU","VAL","ILE")

### --- 3A. FUNZIONE DI ESTRAZIONE GEOMETRIA (eseguita 1 volta per PDB) --- ###

extract_geometry_contacts <- function(pdb_path) {
  tryCatch({
    file_name <- basename(pdb_path)
    
    pdb_aus <- read.pdb(pdb_path)  
    parts <- strsplit(file_name, "_")[[1]]
    chain_HL <- c(parts[2], parts[3])
    chain_AG <- parts[4]
    
    renumbered_df <- renumber_ab_chains(pdb_aus, pdb_path, log_file = file.path(results_dir, "scoring_errors.log"))
    if (!renumbered_df$ok) stop(paste("Renumbering failed:", renumbered_df$error))
    
    df_coord_renumbered <- renumbered_df$df_coord_renumbered
    
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
    
    if (length(BS_A) == 0 || length(BS_HL) == 0) stop("Interfaccia vuota (BS_A o BS_HL è 0)")
    
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
    
    return(contacts) # Restituisce solo il dataframe dei contatti
    
  }, error = function(e) {
    cat("!!! ERRORE su file:", basename(pdb_path), "| Msg:", e$message, "\n")
    return(NULL) # In caso di errore, restituisce NULL
  })
}

### --- 3B. FUNZIONE DI SCORING (eseguita N volte sui contatti già calcolati) --- ###

apply_potential_score <- function(contacts, V_potential_df, potential_type) {
  
  if (is.null(contacts) || nrow(contacts) == 0) {
    return(list(total_potential = NA, mean_potential = NA))
  }
  
  if (potential_type == "sym") {
    contacts <- contacts %>%
      mutate(
        aa1 = sub("_.*", "", res1),
        aa2 = sub("_.*", "", res2),
        pair_key = paste(pmin(aa1, aa2), pmax(aa1, aa2), sep = "_"),
        pot_col = paste0("V", pmin(ring_res1, ring_res2), "_", pmax(ring_res1, ring_res2))
      )
    
    V_potential_lookup <- V_potential_df %>%
      mutate(pair_key = paste(pmin(aa1, aa2), pmax(aa1, aa2), sep = "_"))
    
    contacts <- contacts %>% left_join(V_potential_lookup, by = "pair_key")
    
  } else { # "asym"
    contacts <- contacts %>%
      mutate(
        aa1 = paste0(sub("_.*", "", res1), "_Ag"), 
        aa2 = paste0(sub("_.*", "", res2), "_Ab"), 
        aa_ag = sub("_.*", "", res1), 
        aa_ab = sub("_.*", "", res2), 
        pair_key = paste(aa_ab, aa_ag, sep = "_"),
        pot_col = paste0("V", pmin(ring_res1, ring_res2), "_", pmax(ring_res1, ring_res2))
      )
    
    V_potential_lookup <- V_potential_df %>%
      mutate(pair_key = paste(aa1, aa2, sep = "_"))
    
    contacts <- contacts %>% left_join(V_potential_lookup, by = "pair_key")
  }
  
  if (nrow(contacts) > 0) {
    contacts$potenziale_totale <- sapply(1:nrow(contacts), function(i) {
      col_name <- contacts$pot_col[i]
      if (col_name %in% names(contacts)) {
        return(contacts[i, col_name])
      } else {
        return(NA) 
      }
    })
  } else {
    contacts$potenziale_totale <- numeric(0)
  }
  
  total_potential <- sum(contacts$potenziale_totale, na.rm = TRUE)
  n_contacts <- nrow(contacts)
  mean_potential <- if (n_contacts > 0) total_potential / n_contacts else NA
  
  return(list(total_potential = total_potential, mean_potential = mean_potential))
}


### --- 4. PRE-CALCOLO GEOMETRIA (IL VERO CAMBIAMENTO) --- ###

all_docked_pdbs <- list.files(pdb_dir, pattern = "*.pdb", recursive = TRUE, full.names = TRUE)
cat("Trovate", length(all_docked_pdbs), "pose PDB da analizzare.\n")

cat("--- FASE A: PRE-CALCOLO GEOMETRIA E CONTATTI PER TUTTE LE POSE ---\n")
# future_map calcola la geometria in parallelo e salva i risultati in una lista in memoria RAM
all_pdb_contacts <- future_map(
  all_docked_pdbs,
  extract_geometry_contacts,
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)
# Assegniamo i nomi ai file per recuperarli dopo
names(all_pdb_contacts) <- basename(all_docked_pdbs)


### --- 5. CARICAMENTO POTENZIALI E SCORING --- ###

potential_rds_paths <- list.files(potentials_dir, 
                                  pattern = "potentials_(sym|asym)_wide\\.rds$", 
                                  recursive = TRUE, 
                                  full.names = TRUE)
cat("Trovati", length(potential_rds_paths), "file di potenziali da processare.\n")

# Ciclo sui potenziali
for (potential_path in potential_rds_paths) {
  
  type_tag <- if (str_detect(basename(potential_path), "asym")) "asym" else "sym"
  
  seed_tag <- str_extract(potential_path, "seed_\\d+") 
  split_tag <- str_extract(potential_path, "split_\\d+")
  group_tag <- str_extract(potential_path, "group_\\d+")
  id_tag <- paste(seed_tag, split_tag, group_tag, sep = "_")
  
  output_filename <- paste0("scores_", id_tag, "_", type_tag, ".csv")
  output_filepath <- file.path(results_dir, output_filename)
  
  if (file.exists(output_filepath)) {
    cat("--- SALTO:", output_filename, "(già calcolato) ---\n")
    next 
  }
  
  cat("--- INIZIO SCORING PER:", output_filename, "---\n")
  V_potential_df <- readRDS(potential_path)
  
  # Usiamo purrr::map_dfr (o future_map_dfr) per applicare il potenziale alla lista di contatti precalcolata.
  # Poiché è solo un'operazione sui dataframe, map_dfr è sufficientemente veloce, 
  # ma manteniamo future_map_dfr per sicurezza sui grandi numeri.
  summary_results <- future_map_dfr(
    names(all_pdb_contacts),
    function(pdb_name) {
      contacts <- all_pdb_contacts[[pdb_name]]
      res <- apply_potential_score(contacts, V_potential_df, type_tag)
      data.frame(
        pdb = pdb_name,
        sum_potential = res$total_potential,
        mean_potential = res$mean_potential
      )
    },
    .progress = FALSE, # Falso per non intasare l'output nel ciclo for
    .options = furrr_options(seed = TRUE)
  )
  
  write.csv(summary_results, file = output_filepath, row.names = FALSE)
  cat("--- SALVATO:", output_filename, "---\n")
} 

cat("--- BLOCCO 6: FINITO! Tutti i", length(potential_rds_paths), "potenziali sono stati processati. ---\n")
