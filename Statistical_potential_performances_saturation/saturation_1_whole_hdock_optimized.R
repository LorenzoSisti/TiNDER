### BLOCCO 6: SCORING WHOLE-INTERFACE (PRE-CALCOLO GEOMETRIA) ###

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

# Imposta il parallelismo per future_map
plan(multisession, workers = parallel::detectCores() - 1)

cat("--- BLOCCO 6: SCRIPT DI SCORING WHOLE-INTERFACE INIZIATO ---\n")
cat("Processori usati:", nbrOfWorkers(), "\n")

### --- 2. DEFINIZIONE DEI PERCORSI (INPUT E OUTPUT) --- ###

# Input 1: La cartella con i potenziali WHOLE INTERFACE salvati come V_sym.csv e V_asym.csv
potentials_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali_marzo/potenziali_statistici_whole_interface/"

# Input 2: La cartella con le 22.000 POSE DI DOCKING da analizzare
pdb_dir <- "/Users/lorenzosisti/Downloads/models/"

# Output: La cartella dove salvare i file CSV di risultati
results_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali_marzo/risultati_scoring_finale_HDOCK_WholeInt/"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Parametri globali
DistCutoff <- 8.5
amino_acids <- c("ARG","LYS","ASN","ASP","GLN","GLU","HIS","PRO","TYR","TRP",
                 "SER","THR","GLY","ALA","MET","CYS","PHE","LEU","VAL","ILE")

### --- 3A. FUNZIONE DI ESTRAZIONE GEOMETRIA E CONTATTI --- ###

extract_geometry_contacts_whole_int <- function(pdb_path) {
  tryCatch({
    file_name <- basename(pdb_path)
    
    pdb_aus <- read.pdb(pdb_path)  
    parts <- strsplit(file_name, "_")[[1]]
    chain_HL <- c(parts[2], parts[3])
    chain_AG <- parts[4]
    
    renumbered_df <- renumber_ab_chains(pdb_aus, pdb_path, log_file = file.path(results_dir, "scoring_errors.log"))
    if (!renumbered_df$ok) stop(paste("Renumbering failed:", renumbered_df$error))
    
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
    
    # Se non ci sono contatti, restituiamo NULL per gestirlo nello step successivo
    if (sum(Inter_DistMat_Bin) == 0) return(NULL) 
    
    contact_indices <- which(Inter_DistMat_Bin == 1, arr.ind = TRUE)
    contacts <- data.frame(
      res1 = rownames(Inter_DistMat_Bin)[contact_indices[, 1]],
      res2 = colnames(Inter_DistMat_Bin)[contact_indices[, 2]]
    )
    
    return(contacts)
    
  }, error = function(e) {
    cat("!!! ERRORE su file:", basename(pdb_path), "| Msg:", e$message, "\n")
    return(NULL)
  })
}

### --- 3B. FUNZIONE DI SCORING SUI CONTATTI PRE-CALCOLATI --- ###

apply_potential_score_whole_int <- function(contacts, V_potential_df, potential_type) {
  
  if (is.null(contacts) || nrow(contacts) == 0) {
    return(list(total_potential = 0, mean_potential = 0))
  }
  
  if (potential_type == "sym") {
    contatti_clean <- contacts %>%
      mutate(
        aa1_tmp = sub("_.*", "", res1),
        aa2_tmp = sub("_.*", "", res2),
        key_pair = paste(pmin(aa1_tmp, aa2_tmp), pmax(aa1_tmp, aa2_tmp), sep="_") 
      ) %>%
      group_by(key_pair) %>%
      summarise(n_contacts = n(), .groups = "drop")
    
    merged_res <- V_potential_df %>%
      mutate(key_pair = paste(pmin(aa1, aa2), pmax(aa1, aa2), sep="_")) %>%
      select(key_pair, value) %>%
      distinct(key_pair, .keep_all = TRUE) %>% 
      right_join(contatti_clean, by = "key_pair") %>%
      mutate(
        potenziale_totale = value * n_contacts
      )
    
  } else { # Asym
    contatti_clean <- contacts %>%
      mutate(
        aa1 = sub("_.*", "", res2), # Anticorpo (res2) -> Righe della matrice
        aa2 = sub("_.*", "", res1)  # Antigene (res1)  -> Colonne della matrice
      ) %>%
      group_by(aa1, aa2) %>%
      summarise(n_contacts = n(), .groups = "drop")
    
    merged_res <- V_potential_df %>%
      right_join(contatti_clean, by = c("aa1", "aa2")) %>%
      mutate(
        potenziale_totale = value * n_contacts
      )
  }
  
  total_potential <- sum(merged_res$potenziale_totale, na.rm = TRUE)
  total_observed_contacts <- sum(contatti_clean$n_contacts)
  
  mean_potential <- if (total_observed_contacts > 0) total_potential / total_observed_contacts else NA
  
  return(list(total_potential = total_potential, mean_potential = mean_potential))
}

### --- 4. FASE A: PRE-CALCOLO GEOMETRIA --- ###

all_docked_pdbs <- list.files(pdb_dir, pattern = "*.pdb", recursive = TRUE, full.names = TRUE)
cat("Trovate", length(all_docked_pdbs), "pose PDB da analizzare.\n")

cat("--- FASE A: ESTRAZIONE CONTATTI (WHOLE INTERFACE) IN CORSO... ---\n")
all_pdb_contacts <- future_map(
  all_docked_pdbs,
  extract_geometry_contacts_whole_int,
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)
names(all_pdb_contacts) <- basename(all_docked_pdbs)

### --- 5. FASE B: SCORING CICLICO --- ###

potential_csv_paths <- list.files(potentials_dir, 
                                  pattern = "V_(sym|asym)\\.csv$", 
                                  recursive = TRUE, 
                                  full.names = TRUE)
cat("Trovati", length(potential_csv_paths), "file di potenziali da processare.\n")

for (potential_path in potential_csv_paths) {
  
  type_tag <- if (str_detect(basename(potential_path), "asym")) "asym" else "sym"
  
  seed_tag <- str_extract(potential_path, "seed_\\d+") 
  split_tag <- str_extract(potential_path, "split_\\d+")
  group_tag <- str_extract(dirname(potential_path), "group_\\d+") 
  
  id_tag <- paste(seed_tag, split_tag, group_tag, sep = "_")
  output_filename <- paste0("scores_", id_tag, "_", type_tag, ".csv")
  output_filepath <- file.path(results_dir, output_filename)
  
  if (file.exists(output_filepath)) {
    cat("--- SALTO:", output_filename, "(già calcolato) ---\n")
    next
  }
  
  cat("--- INIZIO CALCOLO PER:", output_filename, "---\n")
  
  # 1. Leggi la matrice 20x20 CSV
  V_mat <- read.csv(potential_path, row.names = 1)
  
  # 2. Trasformala in dataframe lungo (aa1, aa2, value)
  V_potential_df <- expand.grid(aa1 = rownames(V_mat), aa2 = colnames(V_mat), stringsAsFactors = FALSE)
  V_potential_df$value <- as.vector(as.matrix(V_mat))
  
  # Pulizia suffissi eventuali ("ALA_Ab" ecc.)
  V_potential_df$aa1 <- sub("_.*", "", V_potential_df$aa1)
  V_potential_df$aa2 <- sub("_.*", "", V_potential_df$aa2)
  
  # 3. Applicazione del potenziale alla lista di contatti in memoria
  summary_results <- future_map_dfr(
    names(all_pdb_contacts),
    function(pdb_name) {
      contacts <- all_pdb_contacts[[pdb_name]]
      res <- apply_potential_score_whole_int(contacts, V_potential_df, type_tag)
      data.frame(
        pdb = pdb_name,
        sum_potential = res$total_potential,
        mean_potential = res$mean_potential
      )
    },
    .progress = FALSE, # Falso per non sporcare il log durante il loop
    .options = furrr_options(seed = TRUE)
  )
  
  # 4. Salvataggio
  write.csv(summary_results, file = output_filepath, row.names = FALSE)
  cat("--- SALVATO:", output_filename, "---\n")
} 

cat("--- BLOCCO 6 (Whole Interface): FINITO! Processati", length(potential_csv_paths), "potenziali. ---\n")