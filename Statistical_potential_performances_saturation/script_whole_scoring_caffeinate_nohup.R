### BLOCCO 6 (MODIFICATO): SCORING WHOLE-INTERFACE PER TUTTE LE POSE ###

### --- 1. LIBRERIE E SETUP --- ###
library(bio3d)
library(dplyr)
library(tidyr)
library(future)
library(furrr)
library(purrr)
library(stringr)

# Carica le tue funzioni personalizzate (assicurati che il path sia corretto)
source("/Users/lorenzosisti/Documents/Script_ottimizzati_funzioni/functions.R")

# Imposta il parallelismo per future_map (usa tutti i core tranne 1)
plan(multisession, workers = parallel::detectCores() - 1)

cat("--- BLOCCO 6: SCORING WHOLE-INTERFACE INIZIATO ---\n")
cat("Processori usati:", nbrOfWorkers(), "\n")

### --- 2. DEFINIZIONE DEI PERCORSI (INPUT E OUTPUT) --- ###

# Input 1: La cartella con le LISTE dei potenziali (modifica i path dentro i txt se necessario)
potentials_list_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali/potenziali_selezionati_liste/"

# Input 2: La cartella con le 22.000 POSE DI DOCKING
pdb_dir <- "/Users/lorenzosisti/Downloads/models/"

# Output: La cartella dove salvare i risultati
results_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali/risultati_scoring_finale_HDOCK_WholeInt/"
dir.create(results_dir, showWarnings = FALSE)

# Parametri globali
DistCutoff <- 8.5

### --- 3. FUNZIONE DI SCORING RIFATTORIZZATA (WHOLE INTERFACE) --- ###

calculate_whole_interface_score <- function(pdb_path, V_potential_df, potential_type) {
  
  file_name <- basename(pdb_path)
  
  tryCatch({
    
    # --- A. PREPARAZIONE GEOMETRIA (Comune a Sym e Asym) ---
    pdb_aus <- read.pdb(pdb_path)  
    parts <- strsplit(file_name, "_")[[1]]
    chain_HL <- c(parts[2], parts[3])
    chain_AG <- parts[4]
    
    # Rinumerazione
    renumbered_df <- renumber_ab_chains(pdb_aus, pdb_path, log_file = file.path(results_dir, "scoring_errors.log"))
    if (!renumbered_df$ok) {
      stop(paste("Renumbering failed:", renumbered_df$error))
    }
    df_coord_renumbered <- renumbered_df$df_coord_renumbered
    
    # Calcolo centroidi
    centroidi_df <- as.data.frame(
      df_coord_renumbered %>%
        group_by(chain, resno, resid) %>%
        summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
    )
    
    res_names <- paste(centroidi_df$resid, centroidi_df$resno, centroidi_df$chain, sep = "_")
    df_coord_resid_xyz <- centroidi_df[, c("resid", "resno", "x", "y", "z")] 
    rownames(df_coord_resid_xyz) <- res_names
    
    # Calcolo Distanze e Contatti
    DistMat <- as.matrix(dist(df_coord_resid_xyz[, c("x", "y", "z")]))
    
    condA <- !(centroidi_df$chain %in% chain_HL)
    condHL <- centroidi_df$chain %in% chain_HL
    
    Inter_DistMat <- DistMat[condA, condHL]
    Inter_DistMat_Bin <- ifelse(Inter_DistMat <= DistCutoff, 1, 0)
    
    # Check se l'interfaccia esiste
    if (sum(Inter_DistMat_Bin) == 0) {
      return(list(total_potential = 0, mean_potential = 0))
    }
    
    contact_indices <- which(Inter_DistMat_Bin == 1, arr.ind = TRUE)
    contacts <- data.frame(
      res1 = rownames(Inter_DistMat_Bin)[contact_indices[, 1]],
      res2 = colnames(Inter_DistMat_Bin)[contact_indices[, 2]]
    )
    
    # --- B. CALCOLO POTENZIALE (Branching Sym vs Asym) ---
    
    if (potential_type == "sym") {
      
      # Logica tratta da 'assign_sym_pmf_to_docking'
      contatti_clean <- contacts %>%
        mutate(
          aa1 = sub("_.*", "", res1),
          aa2 = sub("_.*", "", res2),
          # IMPORTANTE: Per i simmetrici, ordiniamo alfabeticamente per matchare 
          # chiavi come "ALA_TRP" indipendentemente dall'ordine nel pdb
          key_pair = paste(pmin(aa1, aa2), pmax(aa1, aa2), sep="_") 
        ) %>%
        # Raggruppiamo per coppia unica simmetrica
        group_by(aa1 = pmin(aa1, aa2), aa2 = pmax(aa1, aa2)) %>%
        summarise(n_contacts = n(), .groups = "drop")
      
      # Join con il dataframe del potenziale (V_potential_df)
      # Assumiamo che V_potential_df abbia colonne: aa1, aa2, value
      merged_res <- V_potential_df %>%
        select(aa1, aa2, value) %>%
        left_join(contatti_clean, by = c("aa1", "aa2")) %>%
        mutate(
          n_contacts = replace_na(n_contacts, 0),
          potenziale_totale = value * n_contacts
        )
      
    } else { # Asym
      
      # Logica tratta da 'assign_asym_pmf_to_docking'
      contatti_clean <- contacts %>%
        mutate(
          # res2 è sulle colonne (chain_HL -> Anticorpo -> Ab)
          # res1 è sulle righe (chain_AG -> Antigene -> Ag)
          # Nota: Verifica sempre che Inter_DistMat sia [Ag, Ab]. Nel codice sopra:
          # Inter_DistMat <- DistMat[condA, condHL] -> Righe=Ag, Colonne=Ab.
          
          aa1 = paste0(sub("_.*", "", res2), "_Ab"), # Anticorpo
          aa2 = paste0(sub("_.*", "", res1), "_Ag")  # Antigene
        ) %>%
        group_by(aa1, aa2) %>%
        summarise(n_contacts = n(), .groups = "drop")
      
      # Join con il dataframe del potenziale
      merged_res <- V_potential_df %>%
        select(aa1, aa2, value) %>%
        left_join(contatti_clean, by = c("aa1", "aa2")) %>%
        mutate(
          n_contacts = replace_na(n_contacts, 0),
          potenziale_totale = value * n_contacts
        )
    }
    
    # --- C. CALCOLO SCORE FINALE ---
    total_potential <- sum(merged_res$potenziale_totale, na.rm = TRUE)
    
    # Calcolo numero totale di contatti osservati per la media
    total_observed_contacts <- sum(contatti_clean$n_contacts)
    
    mean_potential <- if (total_observed_contacts > 0) {
      total_potential / total_observed_contacts
    } else {
      0 # O NA, a seconda della preferenza
    }
    
    return(list(
      total_potential = total_potential, 
      mean_potential = mean_potential
    ))
    
  }, error = function(e) {
    cat("!!! ERRORE su file:", basename(pdb_path), "| Msg:", e$message, "\n")
    return(list(total_potential = NA, mean_potential = NA))
  })
}

### --- 4. CARICAMENTO LISTE E CICLO PRINCIPALE --- ###

all_docked_pdbs <- list.files(pdb_dir, pattern = "*.pdb", recursive = TRUE, full.names = TRUE)
cat("Trovate", length(all_docked_pdbs), "pose PDB da analizzare.\n")

potential_list_files <- list.files(potentials_list_dir, pattern = "\\.txt$", full.names = TRUE)
cat("Trovati", length(potential_list_files), "file di liste di potenziali da processare.\n")

for (list_file_path in potential_list_files) {
  
  potential_rds_paths <- readLines(list_file_path)
  
  for (potential_path in potential_rds_paths) {
    
    # --- A. IDENTIFICAZIONE TIPO (Sym/Asym) e ID ---
    type_tag <- if (str_detect(potential_path, "sym")) "sym" else "asym"
    
    # Estrazione ID come nel vecchio script
    seed_tag <- str_extract(potential_path, "seed_\\d+") 
    split_tag <- str_extract(potential_path, "split_\\d+")
    group_tag <- str_extract(potential_path, "group_\\d+")
    id_tag <- paste(seed_tag, split_tag, group_tag, sep = "_")
    
    # Nome file output (identico formato)
    output_filename <- paste0("scores_", id_tag, "_", type_tag, ".csv")
    output_filepath <- file.path(results_dir, output_filename)
    
    if (file.exists(output_filepath)) {
      cat("--- SALTO:", output_filename, "(già calcolato) ---\n")
      next
    }
    
    cat("--- INIZIO CALCOLO PER:", output_filename, "---\n")
    
    # Carica il potenziale (Whole Interface)
    # Si aspetta un DF con colonne: aa1, aa2, value
    V_potential_df <- readRDS(potential_path)
    
    # --- B. ESECUZIONE PARALLELA ---
    
    summary_results <- future_map_dfr(
      all_docked_pdbs,
      function(pdb_path) {
        
        res <- calculate_whole_interface_score(pdb_path, V_potential_df, type_tag)
        
        data.frame(
          pdb = basename(pdb_path),
          sum_potential = res$total_potential,
          mean_potential = res$mean_potential
        )
      },
      .progress = TRUE,
      .options = furrr_options(seed = TRUE)
    )
    
    # --- C. SALVATAGGIO ---
    write.csv(summary_results, file = output_filepath, row.names = FALSE)
    cat("--- SALVATO:", output_filename, "---\n")
    
  } 
}

cat("--- BLOCCO 6 (Whole Interface): FINITO! ---\n")