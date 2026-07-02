################################################################################
# SCRIPT OTTIMIZZATO PER IL CALCOLO DEI POTENZIALI SPECIFICI PER CDR (SYM E ASYM)
# UTILIZZA DATA.TABLE PER ASSOCIARE I CONTATTI ALLA CORRETTA MATRICE CDR
################################################################################

### 1. CARICAMENTO LIBRERIE
pacman::p_load(bio3d, dplyr, future, furrr, purrr, progressr, data.table, stringr)

# Caricamento funzioni custom
# source("/Users/lorenzosisti/Documents/Script_ottimizzati_funzioni/functions.R")
source("/path/to/your/custom/functions/files/functions.R") 

### 2. CONFIGURAZIONE PARALLELIZZAZIONE
plan(multisession, workers = parallel::detectCores() - 1)
handlers(global = TRUE)
handlers("rstudio")

### 3. PARAMETRI GLOBALI E DIRECTORY
pdb_dir <- "/Users/lorenzosisti/Downloads/models"
results_dir <- "/Users/lorenzosisti/Downloads/potenziali_statistici_30_06_hdock"
dir.create(results_dir, showWarnings = FALSE)

DistCutoff <- 8.5

# Ordine degli aminoacidi
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE")

# Recupero di tutti i file PDB
all_docked_pdbs <- list.files(pdb_dir, pattern = "*.pdb", recursive = TRUE, full.names = TRUE)

### 4. PREPARAZIONE MATRICI DEI POTENZIALI (SOLO CDR)
df_asym_long <- fread("/Users/lorenzosisti/Downloads/potenziali_statistici_whole_29_06_data_table_sippl/whole_int_asym_potential.csv")
df_sym_long  <- fread("/Users/lorenzosisti/Downloads/potenziali_statistici_whole_29_06_data_table_sippl/whole_int_sym_potential.csv")

cdr_parts <- c("h1", "h2", "h3", "l1", "l2", "l3")

# ASIMMETRICO: Filtro solo le CDR e aggiungo i suffissi _Ab e _Ag
asym_potentials_cdr <- df_asym_long[part %in% cdr_parts, .(
  part, # Manteniamo la colonna 'part' per il join
  aa1   = paste0(resid_ab, "_Ab"),
  aa2   = paste0(resid_ag, "_Ag"),
  value = potential
)]

# SIMMETRICO: Filtro solo le CDR e "specchio" le coppie non ordinate
df_sym_cdr <- df_sym_long[part %in% cdr_parts]
sym_potentials_cdr <- rbindlist(list(
  df_sym_cdr[, .(part, aa1 = resid_i, aa2 = resid_j, value = potential)],
  df_sym_cdr[resid_i != resid_j, .(part, aa1 = resid_j, aa2 = resid_i, value = potential)]
))

# ==============================================================================
# FUNZIONE CORE: CALCOLO POTENZIALI PER CDR
# ==============================================================================
score_docking_pose_cdr <- function(pdb_path) {
  file_name <- basename(pdb_path)
  
  tryCatch({
    # 1. Estrazione info catene dal filename (es: name_H_L_A.pdb)
    parts <- strsplit(file_name, "_")[[1]]
    chain_H <- parts[2]
    chain_L <- parts[3]
    chain_HL <- c(chain_H, chain_L)
    chain_AG <- sub("\\.pdb$", "", parts[4]) 
    
    # 2. Lettura e Renumbering
    pdb_aus <- read.pdb(pdb_path)
    renumbered_df <- renumber_ab_chains(pdb_aus, pdb_path, log_file = "errors.log")
    
    if (!renumbered_df$ok) return(NULL)
    
    # 3. Calcolo Centroidi via data.table
    dt_coord <- as.data.table(renumbered_df$df_coord_renumbered)
    dt_centroids <- dt_coord[, .(x = mean(x), y = mean(y), z = mean(z)), by = .(chain, resno, resid)]
    
    # 4. Mappatura CDR tramite fcase() sui residui dell'anticorpo
    dt_ab <- dt_centroids[chain %in% chain_HL]
    dt_ab[, cdr := fcase(
      chain == chain_H & resno %in% 26:32, "h1",
      chain == chain_H & resno %in% 52:56, "h2",
      chain == chain_H & resno %in% 95:102, "h3",
      chain == chain_L & resno %in% 24:34, "l1",
      chain == chain_L & resno %in% 50:56, "l2",
      chain == chain_L & resno %in% 89:97, "l3",
      default = NA_character_
    )]
    
    # Teniamo SOLO i residui dell'anticorpo che fanno parte di una CDR
    dt_ab_cdr <- dt_ab[!is.na(cdr)]
    dt_ag <- dt_centroids[!chain %in% chain_HL]
    
    if (nrow(dt_ab_cdr) == 0 | nrow(dt_ag) == 0) return(NULL)
    
    # 5. Calcolo Contatti (Cartesian Join + Distanza)
    dt_ab_cdr[, .dummy := 1L]
    dt_ag[, .dummy := 1L]
    
    dt_contacts <- dt_ab_cdr[dt_ag, on = ".dummy", allow.cartesian = TRUE][
      , dist := sqrt((x - i.x)^2 + (y - i.y)^2 + (z - i.z)^2)
    ][dist <= DistCutoff & resid %in% amino_acids & i.resid %in% amino_acids]
    
    if (nrow(dt_contacts) == 0) {
      message("⚠️ Nessun contatto CDR trovato in ", file_name)
      return(NULL)
    }
    
    # 6. Conteggio coppie specifiche per singola CDR (raggruppo anche per cdr)
    dt_pairs <- dt_contacts[, .(n_contacts = .N), by = .(cdr, resid_ab = resid, resid_ag = i.resid)]
    
    # ---------------------------------------------------------
    # JOIN ASIMMETRICO CON I POTENZIALI DELLE CDR
    # ---------------------------------------------------------
    dt_asym <- copy(dt_pairs)
    dt_asym[, `:=`(part = cdr, aa1 = paste0(resid_ab, "_Ab"), aa2 = paste0(resid_ag, "_Ag"))]
    
    # Faccio il merge considerando "part" (la CDR), "aa1" e "aa2"
    dt_asym <- merge(dt_asym, asym_potentials_cdr, by = c("part", "aa1", "aa2"), all.x = TRUE)
    score_asym <- sum(dt_asym$value * dt_asym$n_contacts, na.rm = TRUE)
    
    # ---------------------------------------------------------
    # JOIN SIMMETRICO CON I POTENZIALI DELLE CDR
    # ---------------------------------------------------------
    dt_sym <- copy(dt_pairs)
    dt_sym[, `:=`(part = cdr, aa1 = resid_ab, aa2 = resid_ag)]
    
    dt_sym <- merge(dt_sym, sym_potentials_cdr, by = c("part", "aa1", "aa2"), all.x = TRUE)
    score_sym <- sum(dt_sym$value * dt_sym$n_contacts, na.rm = TRUE)
    
    # ---------------------------------------------------------
    # OUTPUT
    # ---------------------------------------------------------
    total_contacts <- sum(dt_pairs$n_contacts)
    
    return(data.table(
      pdb_filename      = file_name,
      total_contacts    = total_contacts,
      score_global_sym  = score_sym,
      score_global_asym = score_asym,
      mean_sym          = if(total_contacts > 0) score_sym / total_contacts else NA,
      mean_asym         = if(total_contacts > 0) score_asym / total_contacts else NA
    ))
    
  }, error = function(e) { 
    message("Errore in ", file_name, ": ", e$message)
    return(NULL) 
  })
}

# ==============================================================================
# 5. ESECUZIONE PARALLELA (FUTURE_MAP)
# ==============================================================================
with_progress({
  summary_results_list <- future_map(
    all_docked_pdbs,
    score_docking_pose_cdr,
    .progress = TRUE,
    .options  = furrr_options(seed = TRUE)
  )
})

summary_results <- rbindlist(compact(summary_results_list), use.names = TRUE, fill = TRUE)

### 6. SALVATAGGIO DEI RISULTATI FINALI
fwrite(summary_results, file = file.path(results_dir, "punteggi_cdr_per_posa.csv"))

cat("\n✅ Calcolo dei potenziali CDR completato con successo!\n")
