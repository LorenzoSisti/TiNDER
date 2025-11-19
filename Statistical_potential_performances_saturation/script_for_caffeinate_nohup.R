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
potentials_list_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali/potenziali_selezionati_liste/"

# Input 2: La cartella con le 22.000 POSE DI DOCKING da analizzare
pdb_dir <- "/Users/lorenzosisti/Downloads/models/"

# Output: La cartella dove salvare i 40 file CSV di risultati
results_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali/risultati_scoring_finale_HDOCK/"
dir.create(results_dir, showWarnings = FALSE)

# Parametri globali
DistCutoff <- 8.5
Nsteps <- 3
amino_acids <- c("ARG","LYS","ASN","ASP","GLN","GLU","HIS","PRO","TYR","TRP",
                 "SER","THR","GLY","ALA","MET","CYS","PHE","LEU","VAL","ILE")

# --- 3. FUNZIONE DI SCORING RIFATTORIZZATA --- ###
# Ho unito le tue due funzioni 'assign_sym' e 'assign_asym' in una sola
# per evitare di ripetere il 95% del codice (parsing PDB, centroidi, anelli).

calculate_pose_score <- function(pdb_path, V_potential_df, potential_type) {
  
  # Usiamo un tryCatch per ogni singolo file PDB.
  # Se un PDB è corrotto, logga l'errore e restituisce NA, senza bloccare
  # l'intero calcolo delle altre 21.999 pose.
  tryCatch({
    
    file_name <- basename(pdb_path)
    
    pdb_aus <- read.pdb(pdb_path)  
    parts <- strsplit(file_name, "_")[[1]]
    chain_HL <- c(parts[2], parts[3])
    chain_AG <- parts[4]
    
    # Esegui rinumerazione (dalla tua funzione)
    renumbered_df <- renumber_ab_chains(pdb_aus, pdb_path, log_file = file.path(results_dir, "scoring_errors.log"))
    if (!renumbered_df$ok) {
      stop(paste("Renumbering failed:", renumbered_df$error))
    }
    df_coord_renumbered <- renumbered_df$df_coord_renumbered
    
    # Calcolo centroidi (dalla tua funzione)
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
    
    if (length(BS_A) == 0 || length(BS_HL) == 0) {
      stop("Interfaccia vuota (BS_A o BS_HL è 0)")
    }
    
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
    
    # --- Calcolo potenziale (Logica IF/ELSE) ---
    if (potential_type == "sym") {
      
      contacts <- contacts %>%
        mutate(
          aa1 = sub("_.*", "", res1),
          aa2 = sub("_.*", "", res2),
          # Crea la coppia aa1-aa2 e la colonna potenziale da cercare
          pair_key = paste(pmin(aa1, aa2), pmax(aa1, aa2), sep = "_"),
          pot_col = paste0("V", pmin(ring_res1, ring_res2), "_", pmax(ring_res1, ring_res2))
        )
      
      # Uniamo V_potential_df per efficienza (molto più veloce di mapply)
      V_potential_lookup <- V_potential_df %>%
        mutate(pair_key = paste(pmin(aa1, aa2), pmax(aa1, aa2), sep = "_"))
      
      contacts <- contacts %>%
        left_join(V_potential_lookup, by = "pair_key")
      
    } else { # "asym"
      
      contacts <- contacts %>%
        mutate(
          aa1 = paste0(sub("_.*", "", res1), "_Ag"), # Antigene (res1)
          aa2 = paste0(sub("_.*", "", res2), "_Ab"), # Anticorpo (res2)
          # NOTA: Assumendo res1 = Ag, res2 = Ab. Inverti se necessario.
          # La tua funzione asym usava res1->Ag, res2->Ab. Ho corretto la mia logica
          # in base al tuo codice originale (aa1 = Ab, aa2 = Ag)
          
          # Dalla tua funzione asym:
          aa_ag = sub("_.*", "", res1), # res1 = antigene
          aa_ab = sub("_.*", "", res2), # res2 = anticorpo
          
          pair_key = paste(aa_ab, aa_ag, sep = "_"), # Chiave = Ab_Ag
          pot_col = paste0("V", pmin(ring_res1, ring_res2), "_", pmax(ring_res1, ring_res2))
        )
      
      # Uniamo V_potential_df (che ha colonne aa1=Ab, aa2=Ag)
      V_potential_lookup <- V_potential_df %>%
        mutate(pair_key = paste(aa1, aa2, sep = "_"))
      
      contacts <- contacts %>%
        left_join(V_potential_lookup, by = "pair_key")
    }
    
    # Estrai il valore del potenziale dalla colonna corretta
    # Questo approccio "row/column indexing" è veloce
    if (nrow(contacts) > 0) {
      contacts$potenziale_totale <- sapply(1:nrow(contacts), function(i) {
        col_name <- contacts$pot_col[i]
        if (col_name %in% names(contacts)) {
          return(contacts[i, col_name])
        } else {
          return(NA) # Colonna non trovata (es. V4_4 se Nsteps=3)
        }
      })
    } else {
      contacts$potenziale_totale <- numeric(0)
    }
    
    total_potential <- sum(contacts$potenziale_totale, na.rm = TRUE)
    n_contacts <- nrow(contacts)
    mean_potential <- if (n_contacts > 0) total_potential / n_contacts else NA
    
    return(list(
      total_potential = total_potential,
      mean_potential = mean_potential
    ))
    
  }, error = function(e) {
    # Se questo PDB fallisce, logga e restituisci NA
    cat("!!! ERRORE su file:", basename(pdb_path), "| Msg:", e$message, "\n")
    return(list(total_potential = NA, mean_potential = NA))
  })
}


### --- 4. CARICAMENTO LISTE E CICLO PRINCIPALE --- ###

# Input 1: Lista di tutte le 22.000 pose PDB
all_docked_pdbs <- list.files(pdb_dir, pattern = "*.pdb", recursive = TRUE, full.names = TRUE)
cat("Trovate", length(all_docked_pdbs), "pose PDB da analizzare.\n")

# Input 2: Lista dei 20 file .txt che puntano ai 40 potenziali
potential_list_files <- list.files(potentials_list_dir, pattern = "\\.txt$", full.names = TRUE)
cat("Trovati", length(potential_list_files), "file di liste di potenziali da processare.\n")

# Inizia il ciclo esterno (sui 20 file .txt)
for (list_file_path in potential_list_files) {
  
  potential_rds_paths <- readLines(list_file_path)
  
  # Inizia il ciclo interno (sui 2 percorsi .rds per file .txt)
  for (potential_path in potential_rds_paths) {
    
    # --- A. PREPARAZIONE ---
    
    # Controlla se il file di output esiste GIA'. Se sì, SALTA.
    # Questo rende lo script "resumable"
    type_tag <- if (str_detect(potential_path, "sym")) "sym" else "asym"
    
    # Crea un ID unico per il file di output
    # es: seed_101_split_2_group_1_sym
    # --- NUOVO BLOCCO DI ESTRAZIONE ID ---
    # 1. Estrai il numero di seed
    seed_tag <- str_extract(potential_path, "seed_\\d+") 
    
    # 2. Estrai il numero di split
    split_tag <- str_extract(potential_path, "split_\\d+")
    
    # 3. Estrai il numero di group
    group_tag <- str_extract(potential_path, "group_\\d+")
    
    # 4. Unisci i pezzi in un unico ID (es: seed_123_split_3_group_1)
    # Uso 'paste' per unire i frammenti estratti
    id_tag <- paste(seed_tag, split_tag, group_tag, sep = "_")
    # --- FINE NUOVO BLOCCO DI ESTRAZIONE ID ---
    
    # Il resto del tuo codice per il filename:
    output_filename <- paste0("scores_", id_tag, "_", type_tag, ".csv")
    output_filepath <- file.path(results_dir, output_filename)
    
    if (file.exists(output_filepath)) {
      cat("--- SALTO:", output_filename, "(già calcolato) ---\n")
      next # Salta al prossimo potenziale
    }
    
    cat("--- INIZIO CALCOLO PER:", output_filename, "---\n")
    
    # Carica il potenziale (es. V_sym_wide per questo ciclo)
    V_potential_df <- readRDS(potential_path)
    
    # --- B. ESECUZIONE PARALLELA (su 22.000 PDB) ---
    
    summary_results <- future_map_dfr(
      all_docked_pdbs,
      function(pdb_path) {
        
        # Esegui la funzione di scoring
        res <- calculate_pose_score(pdb_path, V_potential_df, type_tag)
        
        # Restituisci un dataframe riga
        data.frame(
          pdb = basename(pdb_path),
          sum_potential = res$total_potential,
          mean_potential = res$mean_potential
        )
      },
      .progress = TRUE, # Mostra la barra di progresso
      .options = furrr_options(seed = TRUE)
    )
    
    # --- C. SALVATAGGIO ---
    
    write.csv(summary_results, file = output_filepath, row.names = FALSE)
    cat("--- SALVATO:", output_filename, "---\n")
    
  } # Fine ciclo interno (percorsi .rds)
} # Fine ciclo esterno (file .txt)

cat("--- BLOCCO 6: FINITO! Tutti i 40 potenziali sono stati processati. ---\n")