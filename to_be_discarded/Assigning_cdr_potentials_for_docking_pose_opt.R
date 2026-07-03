################################################################################
# SCRIPT OTTIMIZZATO PER IL CALCOLO DEI POTENZIALI STATISTICI (SYM E ASYM)
# STRUTTURATO A FUNZIONI PARALLELIZZATE CON LOGICA CDR E SUFFISSI AB/AG
################################################################################

### 1. CARICAMENTO LIBRERIE
library(bio3d)
library(dplyr)
library(tidyr)
library(readr)
library(future)
library(furrr)

### 2. CONFIGURAZIONE PARALLELIZZAZIONE (Dal tuo secondo script)
plan(multisession, workers = parallel::detectCores() - 1)

### 3. PARAMETRI GLOBALI E DIRECTORY
pdb_dir <- "/Users/lorenzosisti/Downloads/docked_structures_renamed_AF3_11_06"
results_dir <- "/Users/lorenzosisti/Downloads/potenziali_statistici_12_06_50_pose"
dir.create(results_dir, showWarnings = FALSE)

DistCutoff <- 8.5
set.seed(1234)

# Ordine degli aminoacidi (Kyte-Doolittle)
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE") 

# Definizione strutturata delle CDR
cdr_index <- c("cdr_h1", "cdr_h2", "cdr_h3", "cdr_l1", "cdr_l2", "cdr_l3")
cdr_definitions <- list(
  cdr_h1 = list(chain = "H", resno = 26:32),
  cdr_h2 = list(chain = "H", resno = 52:56),
  cdr_h3 = list(chain = "H", resno = 95:102),
  cdr_l1 = list(chain = "L", resno = 24:34),
  cdr_l2 = list(chain = "L", resno = 50:56),
  cdr_l3 = list(chain = "L", resno = 89:97)
)

### 4. CARICAMENTO E PREPARAZIONE MATRICI DEI POTENZIALI (Formato Long)
sym_cdr_potentials  <- read_csv("/Users/lorenzosisti/Downloads/potenziali_statistici_cdr_12_06/V_sym_wide.csv", show_col_types = FALSE)
asym_cdr_potentials <- read_csv("/Users/lorenzosisti/Downloads/potenziali_statistici_cdr_12_06/V_asym_wide.csv", show_col_types = FALSE)

# Trasformazione in formato Long per join rapidi
sym_long <- sym_cdr_potentials %>%
  pivot_longer(cols = starts_with("V_"), names_to = "cdr", names_prefix = "V_", values_to = "potential_sym")

asym_long <- asym_cdr_potentials %>%
  pivot_longer(cols = starts_with("V_"), names_to = "cdr", names_prefix = "V_", values_to = "potential_asym")

# Recupero di tutti i file PDB
pdb_files <- list.files(path = pdb_dir, pattern = "\\.pdb$", full.names = TRUE)


# ==============================================================================
# FUNZIONE 1: CALCOLO POTENZIALI SIMMETRICI PER CDR
# ==============================================================================
assign_sym_cdr_to_docking <- function(pdb_path) {
  
  tryCatch({
    pdb_aus <- read.pdb(pdb_path) 
    
    # Rappresentazione Coarse-grained (Filtro dal tuo primo codice)
    df_coord <- pdb_aus$atom[pdb_aus$atom$type != "HETATM" & (
      (pdb_aus$atom$resid != "GLY" & !(pdb_aus$atom$elety %in% c("N", "CA", "C", "O"))) |
        (pdb_aus$atom$resid == "GLY" & pdb_aus$atom$elety == "CA")),]
    
    if (nrow(df_coord) == 0) return(NULL)
    
    df_coord$insert[is.na(df_coord$insert)] <- ""
    
    # Calcolo coordinate dei centroidi
    centroidi_df <- as.data.frame(
      df_coord %>%
        group_by(chain, resno, insert, resid) %>%
        summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
    )
    
    res_names <- paste(centroidi_df$resid, paste0(centroidi_df$resno, centroidi_df$insert), centroidi_df$chain, sep = "_")
    rownames(centroidi_df) <- res_names
    
    # Matrice di distanza binaria
    DistMat <- as.matrix(dist(centroidi_df[, c("x", "y", "z")])) 
    DistMat_bin <- ifelse(DistMat <= DistCutoff, 1, 0) 
    
    condA <- !(centroidi_df$chain %in% c("H", "L"))
    condHL <- centroidi_df$chain %in% c("H", "L")   
    Inter_DistMat_Bin <- DistMat_bin[condA, condHL, drop = FALSE]
    
    current_pdb_contacts <- list() 
    contact_counter <- 1
    
    # Iterazione sulle 6 CDR
    for (cdr_name in cdr_index) {
      cdr_chain <- cdr_definitions[[cdr_name]]$chain
      cdr_resnos <- cdr_definitions[[cdr_name]]$resno
      
      cond_this_cdr <- centroidi_df$chain == cdr_chain & centroidi_df$resno %in% cdr_resnos
      this_cdr_residue_names <- rownames(centroidi_df)[cond_this_cdr]
      
      valid_cols <- intersect(colnames(Inter_DistMat_Bin), this_cdr_residue_names)
      if (length(valid_cols) == 0) next 
      
      Inter_DistMat_Bin_CDR <- Inter_DistMat_Bin[, valid_cols, drop = FALSE]
      BS_HL_cdr <- apply(Inter_DistMat_Bin_CDR, 2, sum)
      BS_HL_cdr <- BS_HL_cdr[BS_HL_cdr != 0]
      BS_HL_names_cdr <- names(BS_HL_cdr)
      
      for (res_ab in BS_HL_names_cdr) {
        idx_ab <- which(colnames(Inter_DistMat_Bin_CDR) == res_ab)
        contacts <- which(Inter_DistMat_Bin_CDR[, idx_ab] == 1)
        
        for (idx_ag in contacts) {
          res_ag <- rownames(Inter_DistMat_Bin_CDR)[idx_ag]
          aa_ab <- strsplit(res_ab, "_")[[1]][1]
          aa_ag <- strsplit(res_ag, "_")[[1]][1]
          
          if (aa_ab %in% amino_acids && aa_ag %in% amino_acids) {
            current_pdb_contacts[[contact_counter]] <- data.frame(
              cdr = cdr_name, aa1 = aa_ab, aa2 = aa_ag, stringsAsFactors = FALSE
            )
            contact_counter <- contact_counter + 1
          }
        }
      }
    }
    
    if (length(current_pdb_contacts) == 0) return(list(total_potential = 0, n_contacts = 0, mean_potential = 0))
    
    # Collasso e calcolo del potenziale Simmetrico (Nomi puliti, es. "ARG")
    df_scored <- bind_rows(current_pdb_contacts) %>%
      group_by(cdr, aa1, aa2) %>%
      summarise(count = n(), .groups = "drop") %>%
      left_join(sym_long, by = c("cdr", "aa1", "aa2")) %>%
      mutate(total_potential_sym = count * potential_sym)
    
    total_potential <- sum(df_scored$total_potential_sym, na.rm = TRUE)
    n_contacts <- sum(df_scored$count, na.rm = TRUE)
    mean_potential <- if (n_contacts > 0) total_potential / n_contacts else NA
    
    return(list(total_potential = total_potential, n_contacts = n_contacts, mean_potential = mean_potential))
    
  }, error = function(e) { return(NULL) })
}


# ==============================================================================
# FUNZIONE 2: CALCOLO POTENZIALI ASIMMETRICI PER CDR
# ==============================================================================
assign_asym_cdr_to_docking <- function(pdb_path) {
  
  tryCatch({
    pdb_aus <- read.pdb(pdb_path) 
    
    # Rappresentazione Coarse-grained
    df_coord <- pdb_aus$atom[pdb_aus$atom$type != "HETATM" & (
      (pdb_aus$atom$resid != "GLY" & !(pdb_aus$atom$elety %in% c("N", "CA", "C", "O"))) |
        (pdb_aus$atom$resid == "GLY" & pdb_aus$atom$elety == "CA")),]
    
    if (nrow(df_coord) == 0) return(NULL)
    
    df_coord$insert[is.na(df_coord$insert)] <- ""
    
    # Calcolo coordinate dei centroidi
    centroidi_df <- as.data.frame(
      df_coord %>%
        group_by(chain, resno, insert, resid) %>%
        summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
    )
    
    res_names <- paste(centroidi_df$resid, paste0(centroidi_df$resno, centroidi_df$insert), centroidi_df$chain, sep = "_")
    rownames(centroidi_df) <- res_names
    
    # Matrice di distanza binaria
    DistMat <- as.matrix(dist(centroidi_df[, c("x", "y", "z")])) 
    DistMat_bin <- ifelse(DistMat <= DistCutoff, 1, 0) 
    
    condA <- !(centroidi_df$chain %in% c("H", "L"))
    condHL <- centroidi_df$chain %in% c("H", "L")   
    Inter_DistMat_Bin <- DistMat_bin[condA, condHL, drop = FALSE]
    
    current_pdb_contacts <- list() 
    contact_counter <- 1
    
    # Iterazione sulle 6 CDR
    for (cdr_name in cdr_index) {
      cdr_chain <- cdr_definitions[[cdr_name]]$chain
      cdr_resnos <- cdr_definitions[[cdr_name]]$resno
      
      cond_this_cdr <- centroidi_df$chain == cdr_chain & centroidi_df$resno %in% cdr_resnos
      this_cdr_residue_names <- rownames(centroidi_df)[cond_this_cdr]
      
      valid_cols <- intersect(colnames(Inter_DistMat_Bin), this_cdr_residue_names)
      if (length(valid_cols) == 0) next 
      
      Inter_DistMat_Bin_CDR <- Inter_DistMat_Bin[, valid_cols, drop = FALSE]
      BS_HL_cdr <- apply(Inter_DistMat_Bin_CDR, 2, sum)
      BS_HL_cdr <- BS_HL_cdr[BS_HL_cdr != 0]
      BS_HL_names_cdr <- names(BS_HL_cdr)
      
      for (res_ab in BS_HL_names_cdr) {
        idx_ab <- which(colnames(Inter_DistMat_Bin_CDR) == res_ab)
        contacts <- which(Inter_DistMat_Bin_CDR[, idx_ab] == 1)
        
        for (idx_ag in contacts) {
          res_ag <- rownames(Inter_DistMat_Bin_CDR)[idx_ag]
          aa_ab <- strsplit(res_ab, "_")[[1]][1]
          aa_ag <- strsplit(res_ag, "_")[[1]][1]
          
          if (aa_ab %in% amino_acids && aa_ag %in% amino_acids) {
            current_pdb_contacts[[contact_counter]] <- data.frame(
              cdr = cdr_name, aa1 = aa_ab, aa2 = aa_ag, stringsAsFactors = FALSE
            )
            contact_counter <- contact_counter + 1
          }
        }
      }
    }
    
    if (length(current_pdb_contacts) == 0) return(list(total_potential = 0, n_contacts = 0, mean_potential = 0))
    
    # --- LA TUA IDENTICA LOGICA CON SUFFISSI ---
    # Sovrascriviamo aa1 e aa2 incollando subito i suffissi prima del join asimmetrico
    df_scored <- bind_rows(current_pdb_contacts) %>%
      group_by(cdr, aa1, aa2) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(
        aa1 = paste0(aa1, "_Ab"),  # Aggiunge suffisso Anticorpo
        aa2 = paste0(aa2, "_Ag")   # Aggiunge suffisso Antigene
      ) %>%
      left_join(asym_long, by = c("cdr", "aa1", "aa2")) %>%
      mutate(total_potential_asym = count * potential_asym)
    
    total_potential <- sum(df_scored$total_potential_asym, na.rm = TRUE)
    n_contacts <- sum(df_scored$count, na.rm = TRUE)
    mean_potential <- if (n_contacts > 0) total_potential / n_contacts else NA
    
    return(list(total_potential = total_potential, n_contacts = n_contacts, mean_potential = mean_potential))
    
  }, error = function(e) { return(NULL) })
}


# ==============================================================================
# 5. ESECUZIONE PARALLELA (FUTURE_MAP_DFR) - IDENTICO AL SECONDO CODICE
# ==============================================================================
summary_results <- future_map_dfr(
  pdb_files,
  function(pdb_path) {
    message("Elaborazione file: ", basename(pdb_path))
    
    # Chiamata alle due nuove funzioni indipendenti
    res_sym  <- assign_sym_cdr_to_docking(pdb_path)
    res_asym <- assign_asym_cdr_to_docking(pdb_path)
    
    # Gestione dei casi di errore o file vuoti
    if (is.null(res_sym) | is.null(res_asym)) {
      return(data.frame(
        pdb_filename    = basename(pdb_path),
        total_contacts  = NA,
        score_global_sym  = NA,
        score_global_asym = NA,
        mean_sym        = NA,
        mean_asym       = NA
      ))
    }
    
    # Generazione della riga di output per il PDB corrente
    data.frame(
      pdb_filename      = basename(pdb_path),
      total_contacts    = res_sym$n_contacts,
      score_global_sym  = res_sym$total_potential,
      score_global_asym = res_asym$total_potential,
      mean_sym          = res_sym$mean_potential,
      mean_asym         = res_asym$mean_potential
    )
  },
  .progress = TRUE,
  .options  = furrr_options(seed = TRUE)
)

### 6. SALVATAGGIO DEI RISULTATI FINALI
write_csv(summary_results, file = file.path(results_dir, "punteggi_cdr_per_posa.csv"))

cat("Calcolo parallelo a funzioni completato con successo!\n")
