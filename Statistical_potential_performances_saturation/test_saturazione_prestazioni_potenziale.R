### Required libraries
library(bio3d)
library(dplyr)
library(future)
library(furrr)
library(purrr)
library(pheatmap)
library(ggplotify)
library(reshape2)
library(tidyr)
library(tidyverse)

# Define the path to a custom function files
source("/Users/lorenzosisti/Documents/Script_ottimizzati_funzioni/functions.R")

### Set up parallelization to speed up computation
plan(multisession, workers = parallel::detectCores() - 1)

### Define directories and global parameters
pdb_dir <- "/Users/lorenzosisti/Downloads/database_settembre_renamed/"

DistCutoff <- 8.5  
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE")

results_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali"
dir.create(results_dir, showWarnings = FALSE)

seed_list <- c(101, 123, 202, 303, 404, 505, 606, 707, 808, 909)
partitioning_numbers <- c(2, 3, 5, 8, 13, 21, 34, 55, 89, 144)
pdbL <- list.files(pdb_dir, pattern = "\\.pdb$", full.names = TRUE)


for (i in seq_along(seed_list)) {
  
  current_seed <- seed_list[i]
  
  for (n in partitioning_numbers) {
    
    set.seed(current_seed)
    shuffled_files <- sample(pdbL)
    
    # Calcola quanti file per gruppo (floor)
    files_per_group <- floor(length(shuffled_files) / n)
    
    # Prendi solo i file che useremo (gli altri verranno ignorati)
    usable_files <- shuffled_files[1:(files_per_group * n)]
    
    # Crea i gruppi (questa è una lista di vettori di percorsi)
    split_groups <- split(usable_files, rep(1:n, each = files_per_group))
    
    # Ogni run con seed diversa avrà una directory dedicata
    # --- CORREZIONE: Usiamo 'results_dir' invece di 'output_base_dir' ---
    split_dir <- file.path(results_dir, paste0("seed_", current_seed), paste0("split_", n))
    dir.create(split_dir, recursive = TRUE, showWarnings = FALSE)
    
    #
    # --- ### INIZIO MODIFICA PRINCIPALE ### ---
    #
    # Invece di creare una cartella "group_j" e copiarci i file,
    # creiamo un file "group_j.txt" e ci scriviamo dentro la lista dei percorsi.
    #
    for (j in seq_along(split_groups)) {
      
      # 1. Definisci il NOME del file di output (es. group_1.txt)
      output_list_file <- file.path(split_dir, paste0("group_", j, ".txt"))
      
      # 2. Estrai l'elenco dei file PDB per questo gruppo
      files_for_this_group <- split_groups[[j]]
      
      # 3. Scrivi quell'elenco nel file di testo, un percorso per riga
      writeLines(files_for_this_group, con = output_list_file)
    }
    #
    # --- ### FINE MODIFICA PRINCIPALE ### ---
    #
    
    # Salva il seed usato per tracciabilità
    cat("Seed usata:", current_seed, "\n", file = file.path(split_dir, "seed_used.txt"))
    
    # Messaggio di stato aggiornato
    cat("Create", n, "liste di file nella cartella:", split_dir, "\n")
    
  }
}

################################################################################
### SECONDO BLOCCO DI CODICE: CALCOLARE LE MATRICI DI CONTATTO WHOLE-INTERFACE ###
################################################################################

# --- 1. Definiamo chiaramente le cartelle di IN e OUT ---

# Questa è la cartella DOVE IL BLOCCO 1 HA SALVATO LE LISTE .txt
list_base_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali"

# Questa è la cartella DOVE VOGLIAMO SALVARE LE MATRICI FINALI
matrix_output_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali/contact_matrices_dir/"
dir.create(matrix_output_dir, showWarnings = FALSE)

set.seed(1234) # Per la riproducibilità

# --- 2. Troviamo i file .txt, non le cartelle ---
list_files <- list.files(list_base_dir, 
                         pattern = "^group_.*\\.txt$", 
                         recursive = TRUE, 
                         full.names = TRUE)

cat("Trovati", length(list_files), "file di gruppi (.txt) da processare.\n")


# Definiamo una funzione generale per elaborare ogni file PDB
compute_whole_int_statistical_potential_saturation <- function(pdb_path) {
  tryCatch({
    
    # --- Inizio sezione invariata (lettura, rinumerazione, centroidi) ---
    pdb_aus <- read.pdb(pdb_path)
    chain_HL <- c("H", "L")
    chain_AG <- "A"
    
    renumbered_df <- renumber_ab_chains(pdb_aus, pdb_path, log_file = "errors.log")
    if (!renumbered_df$ok) {
      # --- MODIFICA: Miglior log di errore ---
      return(list(ok = FALSE, error = paste("Errore in renumber_ab_chains:", renumbered_df$error, "File:", pdb_path)))
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
    DistMat_bin <- ifelse(DistMat <= DistCutoff, 1, 0)
    
    condA <- !(centroidi_df$chain %in% c("H", "L"))
    condHL <- centroidi_df$chain %in% c("H", "L")
    
    Inter_DistMat <- DistMat[condA, condHL]
    Inter_DistMat_Bin <- ifelse(Inter_DistMat <= DistCutoff, 1, 0)
    
    BS_A <- apply(Inter_DistMat_Bin, 1, sum)
    BS_A <- BS_A[BS_A != 0]
    BS_HL <- apply(Inter_DistMat_Bin, 2, sum)
    BS_HL <- BS_HL[BS_HL != 0]
    
    BS_A_names <- names(BS_A)
    BS_HL_names <- names(BS_HL)
    # --- Fine sezione invariata ---
    
    # --- MODIFICA: Inizializziamo DUE matrici di contatto ---
    contact_matrix_asym <- matrix(0, nrow=20, ncol=20, dimnames=list(amino_acids, amino_acids))
    contact_matrix_sym  <- matrix(0, nrow=20, ncol=20, dimnames=list(amino_acids, amino_acids))
    
    # Inizializziamo i vettori di conteggio (come prima)
    residue_counts_dataset_ab <- setNames(rep(0, length(amino_acids)), amino_acids)
    residue_counts_dataset_ligando <- setNames(rep(0, length(amino_acids)), amino_acids)
    residue_counts_interface_ab <- setNames(rep(0, length(amino_acids)), amino_acids)
    residue_counts_interface_ligando <- setNames(rep(0, length(amino_acids)), amino_acids)
    
    # --- MODIFICA: Riempiamo le matrici (logica presa dal tuo Script 2) ---
    for (res_ab in BS_HL_names) {
      idx_ab <- which(colnames(Inter_DistMat_Bin) == res_ab)
      contacts <- which(Inter_DistMat_Bin[, idx_ab] == 1)
      for (idx_ag in contacts) {
        res_ag <- rownames(Inter_DistMat_Bin)[idx_ag]
        aa_ab <- strsplit(res_ab, "_")[[1]][1]
        aa_ag <- strsplit(res_ag, "_")[[1]][1]
        if (aa_ab %in% amino_acids && aa_ag %in% amino_acids) {
          
          # 1. Matrice Asimmetrica (Righe=Ab, Colonne=Ag)
          contact_matrix_asym[aa_ab, aa_ag] <- contact_matrix_asym[aa_ab, aa_ag] + 2
          
          # 2. Matrice Simmetrica
          if (aa_ab == aa_ag) {
            contact_matrix_sym[aa_ab, aa_ag] <- contact_matrix_sym[aa_ab, aa_ag] + 2
          } else {
            contact_matrix_sym[aa_ab, aa_ag] <- contact_matrix_sym[aa_ab, aa_ag] + 2
            contact_matrix_sym[aa_ag, aa_ab] <- contact_matrix_sym[aa_ag, aa_ab] + 2
          }
        }
      }
    }
    
    # --- Sezione invariata: Riempiamo i vettori di conteggio ---
    # (Questa logica è corretta e calcola sia _dataset_ che _interface_)
    df_centroidi_anticorpo <- centroidi_df[centroidi_df$chain %in% chain_HL, ]
    rownames(df_centroidi_anticorpo) <- paste(df_centroidi_anticorpo$resid, df_centroidi_anticorpo$resno, df_centroidi_anticorpo$chain, sep = "_")
    df_centroidi_anticorpo_bs <- df_centroidi_anticorpo[rownames(df_centroidi_anticorpo) %in% BS_HL_names, ]
    
    residue_counts_interface_ab <- residue_counts_interface_ab + table(factor(df_centroidi_anticorpo_bs$resid, levels = amino_acids))
    residue_counts_dataset_ab <- residue_counts_dataset_ab + table(factor(df_centroidi_anticorpo$resid, levels = amino_acids))
    
    df_centroidi_ligando <- centroidi_df[centroidi_df$chain %in% chain_AG, ]
    rownames(df_centroidi_ligando) <- paste(df_centroidi_ligando$resid, df_centroidi_ligando$resno, df_centroidi_ligando$chain, sep = "_")
    df_centroidi_ligando_bs <- df_centroidi_ligando[rownames(df_centroidi_ligando) %in% BS_A_names, ]
    
    residue_counts_interface_ligando <- residue_counts_interface_ligando + table(factor(df_centroidi_ligando_bs$resid, levels = amino_acids))
    residue_counts_dataset_ligando <- residue_counts_dataset_ligando + table(factor(df_centroidi_ligando$resid, levels = amino_acids))
    
    # --- MODIFICA: Aggiorniamo il return ---
    return(list(
      ok = TRUE,
      contact_matrix_asym = contact_matrix_asym, # Nuovo
      contact_matrix_sym = contact_matrix_sym,   # Sostituisce contact_matrix
      residue_counts_dataset_ab = residue_counts_dataset_ab,
      residue_counts_dataset_ligando = residue_counts_dataset_ligando,
      residue_counts_interface_ab = residue_counts_interface_ab,
      residue_counts_interface_ligando = residue_counts_interface_ligando
    ))
  }, error = function(e) {
    return(list(ok = FALSE, error = paste("Errore generico su file:", pdb_path, "Message:", e$message)))
  })
}

# --- Ciclo principale adattato ---
for (list_file_path in list_files) {
  
  cat("Processando gruppo (da lista):", list_file_path, "\n")
  
  pdb_files <- readLines(list_file_path)
  pdb_files <- pdb_files[pdb_files != ""] 
  
  if (length(pdb_files) == 0) {
    cat("   ... Lista vuota, salto.\n")
    next
  }
  
  results_list <- future_map(pdb_files, compute_whole_int_statistical_potential_saturation, .progress = TRUE)
  
  valid_results <- keep(results_list, ~ .x$ok)
  failed <- keep(results_list, ~ !.x$ok)
  
  # --- Creazione cartella di output ---
  relative_path <- gsub(paste0("^", list_base_dir, "/"), "", list_file_path)
  output_rel_dir <- gsub("\\.txt$", "", relative_path)
  out_dir <- file.path(matrix_output_dir, output_rel_dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Salva i log degli errori
  if (length(failed) > 0) {
    error_messages <- map_chr(failed, "error")
    writeLines(error_messages, file.path(out_dir, "errors.log"))
  }
  
  if (length(valid_results) == 0) {
    cat("   ... Nessun risultato valido per questo gruppo. Creato file di log.\n")
    next
  }
  
  # --- MODIFICA: Aggregazione dei risultati ---
  # Aggreghiamo entrambe le matrici
  contact_matrix_asym <- Reduce(`+`, map(valid_results, "contact_matrix_asym"))
  contact_matrix_sym  <- Reduce(`+`, map(valid_results, "contact_matrix_sym"))
  
  # I conteggi rimangono invariati
  residue_counts_dataset_ab <- Reduce(`+`, map(valid_results, "residue_counts_dataset_ab"))
  residue_counts_dataset_ligando <- Reduce(`+`, map(valid_results, "residue_counts_dataset_ligando"))
  residue_counts_interface_ab <- Reduce(`+`, map(valid_results, "residue_counts_interface_ab"))
  residue_counts_interface_ligando <- Reduce(`+`, map(valid_results, "residue_counts_interface_ligando"))
  
  # --- MODIFICA: Salvataggio di entrambe le matrici ---
  write.csv(contact_matrix_asym, file = file.path(out_dir, "contact_matrix_asym.csv"), row.names = TRUE)
  saveRDS(contact_matrix_asym, file = file.path(out_dir, "contact_matrix_asym.rds"))
  
  write.csv(contact_matrix_sym, file = file.path(out_dir, "contact_matrix_sym.csv"), row.names = TRUE)
  saveRDS(contact_matrix_sym, file = file.path(out_dir, "contact_matrix_sym.rds"))
  
  # Salvataggio dei conteggi (invariato)
  write.csv(residue_counts_dataset_ab, file = file.path(out_dir, "residue_counts_dataset_ab.csv"), row.names = TRUE)
  saveRDS(residue_counts_dataset_ab, file = file.path(out_dir, "residue_counts_dataset_ab.rds"))
  
  write.csv(residue_counts_dataset_ligando, file = file.path(out_dir, "residue_counts_dataset_ligando.csv"), row.names = TRUE)
  saveRDS(residue_counts_dataset_ligando, file = file.path(out_dir, "residue_counts_dataset_ligando.rds"))
  
  write.csv(residue_counts_interface_ab, file = file.path(out_dir, "residue_counts_interface_ab.csv"), row.names = TRUE)
  saveRDS(residue_counts_interface_ab, file = file.path(out_dir, "residue_counts_interface_ab.rds"))
  
  write.csv(residue_counts_interface_ligando, file = file.path(out_dir, "residue_counts_interface_ligando.csv"), row.names = TRUE)
  saveRDS(residue_counts_interface_ligando, file = file.path(out_dir, "residue_counts_interface_ligando.rds"))
  
  cat("   ... Salvataggio completato in:", out_dir, "\n")
}

cat("Finito!\n")

################################################################################
### BLOCCO 2.1: CALCOLARE I POTENZIALI WHOLE-INTERFACE ###
################################################################################

### Directories and global variables
base_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali/contact_matrices_dir/"
output_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali/potenziali_statistici_whole_interface/"
dir.create(output_dir, showWarnings = FALSE)

amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE")

split_dirs <- list.dirs(base_dir, recursive = FALSE)

# Cerca i file "contact_matrix_sym.rds" (salvati dal Blocco 1)
group_dirs <- list.files(base_dir, pattern = "contact_matrix_sym\\.rds", recursive = TRUE, full.names = TRUE) %>%
  dirname() %>%
  unique()

# Lista per raccogliere i risultati emi_matrix
emi_sums <- data.frame(
  group_path = character(),
  emi_sum = numeric(),
  stringsAsFactors = FALSE
)

for (group_path in group_dirs) {
  
  cat("Elaborando:", group_path, "\n")
  
  # Carica matrici
  sym_contact_matrix <- readRDS(file.path(group_path, "contact_matrix_sym.rds"))
  asym_contact_matrix <- readRDS(file.path(group_path, "contact_matrix_asym.rds"))
  residue_counts_interface_ab <- readRDS(file.path(group_path, "residue_counts_interface_ab.rds"))
  residue_counts_interface_ligando <- readRDS(file.path(group_path, "residue_counts_interface_ligando.rds"))
  
  # Vado a definire due quantità per l'approccio asimmetrico
  
  paratope_freq <- residue_counts_interface_ab / sum(residue_counts_interface_ab)
  epitope_freq <- residue_counts_interface_ligando / sum(residue_counts_interface_ligando)
  
  # Ora lo faccio per l'approccio simmetrico
  
  par_plus_epi <- residue_counts_interface_ab + residue_counts_interface_ligando
  residue_freq <- par_plus_epi / sum(par_plus_epi)
  
  # Ora mi calcolo le frequenze di contatto relativo dalla matrice dei contatti
  # Data la simmetria rispetto alla diagonale delle matrici, tutta l'informazione è in una triangolare 
  
  # Compute normalized contact frequencies
  sym_lower_tri <- sym_contact_matrix[lower.tri(sym_contact_matrix, diag = TRUE)] # All the symmetrical contact information is in a triangular
  contact_freq_sym <- sym_contact_matrix / sum(sym_lower_tri) 
  contact_freq_asym <- asym_contact_matrix / sum(asym_contact_matrix)
  rel_path <- gsub(base_dir, "", group_path)
  emi_sums <- rbind(emi_sums, data.frame(group_path = rel_path, emi_sum = sum(sym_lower_tri)))
  
  # Potenziali
  V_asym <- -log(contact_freq_asym / outer(paratope_freq, epitope_freq, "*")) * 2.479
  V_sym <- -log(contact_freq_sym / outer(residue_freq, residue_freq, "*")) * 2.479
  
  V_asym[!is.finite(V_asym)] <- 0
  V_sym[!is.finite(V_sym)] <- 0
  
  # Salva i risultati
  out_subdir <- file.path(output_dir, rel_path)
  dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(V_asym, file = file.path(out_subdir, "V_asym.csv"))
  write.csv(V_sym, file = file.path(out_subdir, "V_sym.csv"))
  
  ### Prepare and plot heatmaps
  # Compute global min/max across both matrices for consistent color scaling
  pot_min <- abs(min(min(V_sym), min(V_asym)))
  pot_max <- abs(max(max(V_sym), max(V_asym)))
  
  # HEATMAP ASIMMETRICA
  rownames(V_asym) <- paste(amino_acids, "Ab", sep = "_")
  colnames(V_asym) <- paste(amino_acids, "Ag", sep = "_")
  
  pdf(file = file.path(out_subdir, "V_asym.pdf"), width = 7, height = 7)
  pheatmap(V_asym, 
           color = colorRampPalette(c("yellow", "green", "cyan", "blue"))(50),  
           breaks = seq(- max(pot_min, pot_max), max(pot_min, pot_max), length.out = 51),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           main = "Asymmetric Pairwise Statistical Potential",
           display_numbers = FALSE,
           fontsize = 10)
  dev.off()
  
  # HEATMAP SIMMETRICA
  pdf(file = file.path(out_subdir, "V_sym.pdf"), width = 7, height = 7)
  pheatmap(V_sym, 
           color = colorRampPalette(c("yellow", "green", "cyan", "blue"))(50),  
           breaks = seq(- max(pot_min, pot_max), max(pot_min, pot_max), length.out = 51),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           main = "Symmetric Pairwise Statistical Potential",
           display_numbers = FALSE,
           fontsize = 10)
  dev.off()
  
  
  cat("Salvato:", out_subdir, "\n")
}

write.csv(emi_sums, file = file.path(output_dir, "emi_matrix_sums.csv"), row.names = FALSE)


################################################################################
### TERZO BLOCCO DI CODICE: CALCOLARE LE MATRICI DI CONTATTO RADIALI ###
################################################################################

# --- 1. Definiamo le cartelle di IN e OUT per questo blocco ---

# Input: Le liste .txt create dal Blocco 1
# (usiamo la STESSA cartella di base dei blocchi precedenti)
list_base_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali"

# Output: Una NUOVA cartella per i risultati radiali (stratificati)
radial_matrix_output_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali/matrici_stratificate/"
dir.create(radial_matrix_output_dir, showWarnings = FALSE)

# --- 2. Parametri specifici per questo blocco ---
Nsteps <- 3 # Numero di "anelli" concentrici
ring_pairs <- c("1-1", "2-2", "3-3", "1-2", "2-3", "1-3")

# --- 3. Troviamo i file .txt da processare ---
# (Questo è lo stesso elenco del Blocco 2)
list_files_to_process <- list.files(list_base_dir, 
                                    pattern = "^group_.*\\.txt$",
                                    recursive = TRUE, 
                                    full.names = TRUE)

cat("\n--- BLOCCO 3 INIZIATO ---\n")
cat("Trovati", length(list_files_to_process), "file di gruppi (.txt) da processare per analisi RADIALE.\n")

# --- 4. Strutture di base per l'inizializzazione ---
base_matrix_20x20 <- matrix(0, nrow = 20, ncol = 20, dimnames = list(amino_acids, amino_acids))
base_vector_20 <- setNames(rep(0, length(amino_acids)), amino_acids)


# --- 5. Funzione di processamento PDB (adattata dal tuo snippet) ---
process_pdb_radial <- function(pdb_path) {
  tryCatch({
    
    # --- I/O e setup (MODIFICATO: usa pdb_path) ---
    pdb_aus <- read.pdb(pdb_path) 
    chain_HL <- c("H", "L")
    chain_AG <- "A"
    
    # --- Coarse-graining (dal tuo snippet) ---
    df_coord <- pdb_aus$atom[pdb_aus$atom$type != "HETATM" & (
      (pdb_aus$atom$resid != "GLY" & !(pdb_aus$atom$elety %in% c("N", "CA", "C", "O"))) |
        (pdb_aus$atom$resid == "GLY" & pdb_aus$atom$elety == "CA")),]
    
    if (nrow(df_coord) == 0) {
      # MODIFICATO: return con 'ok = FALSE'
      return(list(ok = FALSE, error = paste("df_coord vuoto:", pdb_path)))
    }
    
    # --- Renumbering (dal tuo snippet) ---
    df_coord_corrected <- df_coord
    corrected_dfs <- list()
    for (chain in unique(df_coord_corrected$chain)) {
      unique_residues <- unique(paste(df_coord_corrected$resno[df_coord_corrected$chain == chain],
                                      df_coord_corrected$insert[df_coord_corrected$chain == chain], sep=""))
      new_numbering <- setNames(seq_along(unique_residues), unique_residues)
      df_coord_corrected$resno[df_coord_corrected$chain == chain] <- new_numbering[paste(df_coord_corrected$resno[df_coord_corrected$chain == chain],
                                                                                         df_coord_corrected$insert[df_coord_corrected$chain == chain], sep="")]
      corrected_dfs[[chain]] <- df_coord_corrected[df_coord_corrected$chain == chain, ]
    }
    df_coord_corrected <- do.call(rbind, corrected_dfs)
    
    # --- Calcolo Centroidi (dal tuo snippet) ---
    centroids_df <- as.data.frame(
      df_coord_corrected %>%
        group_by(chain, resno, resid) %>%
        summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
    )
    res_names <- paste(centroids_df$resid, centroids_df$resno, centroids_df$chain, sep = "_")
    df_coord_resid_xyz <- centroids_df[, c("resid", "resno", "x", "y", "z")]
    rownames(df_coord_resid_xyz) <- res_names
    DistMat <- as.matrix(dist(df_coord_resid_xyz[, c("x", "y", "z")])) 
    
    # --- Definizione Interfaccia e Anelli (dal tuo snippet) ---
    condA <- !(centroids_df$chain %in% chain_HL)
    condHL <- centroids_df$chain %in% chain_HL
    
    Inter_DistMat <- DistMat[condA, condHL]
    Inter_DistMat_Bin <- ifelse(Inter_DistMat <= DistCutoff, 1, 0)
    
    BS_A <- apply(Inter_DistMat_Bin, 1, sum)
    BS_A <- BS_A[BS_A != 0]
    BS_HL <- apply(Inter_DistMat_Bin, 2, sum)
    BS_HL <- BS_HL[BS_HL != 0]
    
    if (length(BS_A) == 0 || length(BS_HL) == 0) {
      # MODIFICATO: return con 'ok = FALSE'
      return(list(ok = FALSE, error = paste("Interfaccia vuota (BS_A o BS_HL è 0):", pdb_path)))
    }
    
    BS_A_names <- names(BS_A)
    BS_HL_names <- names(BS_HL)
    
    # --- Logica degli anelli (dal tuo snippet) ---
    df_centroids_BS <- df_coord_resid_xyz[rownames(df_coord_resid_xyz) %in% c(BS_A_names, BS_HL_names), ]
    center_BS <- colMeans(df_centroids_BS[, c("x", "y", "z")])
    DistMat_centroid <- as.matrix(dist(rbind(df_centroids_BS[, c("x", "y", "z")], center_BS)))
    vet_dist_centroid <- DistMat_centroid[-nrow(DistMat_centroid), ncol(DistMat_centroid)]
    r_max <- max(vet_dist_centroid)
    borders <- c(0, 0.37, 0.64, 1) * r_max 
    
    # --- Assegnazione contatti e anelli (dal tuo snippet) ---
    contact_indices <- which(Inter_DistMat_Bin == 1, arr.ind = TRUE)
    contacts <- data.frame(
      res1 = rownames(Inter_DistMat_Bin)[contact_indices[, 1]],
      res2 = colnames(Inter_DistMat_Bin)[contact_indices[, 2]]
    )
    contacts$dist_res1_centroid <- vet_dist_centroid[contacts$res1]
    contacts$dist_res2_centroid <- vet_dist_centroid[contacts$res2]
    
    assign_ring <- function(dist) findInterval(dist, borders, left.open = TRUE)
    contacts$ring_res1 <- sapply(contacts$dist_res1_centroid, assign_ring)
    contacts$ring_res2 <- sapply(contacts$dist_res2_centroid, assign_ring)
    
    # --- Inizializzazione strutture LOCALI (per questo PDB) ---
    local_matrices_asym <- setNames(vector("list", length(ring_pairs)), ring_pairs)
    local_matrices_sym <- setNames(vector("list", length(ring_pairs)), ring_pairs)
    for (pair in ring_pairs) {
      local_matrices_asym[[pair]] <- base_matrix_20x20
      local_matrices_sym[[pair]] <- base_matrix_20x20
    }
    
    local_ab_counts <- setNames(vector("list", Nsteps), 1:Nsteps)
    local_ag_counts <- setNames(vector("list", Nsteps), 1:Nsteps)
    for (step in 1:Nsteps) {
      local_ab_counts[[step]] <- base_vector_20
      local_ag_counts[[step]] <- base_vector_20
    }
    
    # --- Riempimento matrici di contatto (dal tuo snippet) ---
    for (k in seq_len(nrow(contacts))) {
      aa1 <- strsplit(contacts$res1[k], "_")[[1]][1]
      aa2 <- strsplit(contacts$res2[k], "_")[[1]][1]
      if (!(aa1 %in% amino_acids && aa2 %in% amino_acids)) next
      
      ring1 <- contacts$ring_res1[k]
      ring2 <- contacts$ring_res2[k]
      
      if (ring1 %in% 1:Nsteps && ring2 %in% 1:Nsteps) {
        rings_sorted <- sort(c(ring1, ring2))
        pair_key <- paste(rings_sorted[1], rings_sorted[2], sep = "-")
        
        chain1 <- strsplit(contacts$res1[k], "_")[[1]][3]
        chain2 <- strsplit(contacts$res2[k], "_")[[1]][3]
        
        # ASYMMETRIC
        if (chain1 %in% chain_HL && chain2 %in% chain_AG) {
          local_matrices_asym[[pair_key]][aa1, aa2] <- local_matrices_asym[[pair_key]][aa1, aa2] + 2
        } else if (chain1 %in% chain_AG && chain2 %in% chain_HL) {
          local_matrices_asym[[pair_key]][aa2, aa1] <- local_matrices_asym[[pair_key]][aa2, aa1] + 2
        }
        
        # SYMMETRIC
        if (aa1 == aa2) {
          local_matrices_sym[[pair_key]][aa1, aa1] <- local_matrices_sym[[pair_key]][aa1, aa1] + 2
        } else {
          local_matrices_sym[[pair_key]][aa1, aa2] <- local_matrices_sym[[pair_key]][aa1, aa2] + 2
          local_matrices_sym[[pair_key]][aa2, aa1] <- local_matrices_sym[[pair_key]][aa2, aa1] + 2
        }
      }
    }
    
    # --- Conteggio residui per anello (dal tuo snippet) ---
    for (step in 1:Nsteps) {
      r_start <- borders[step]
      r_end <- borders[step + 1]
      selected_residues <- names(vet_dist_centroid[vet_dist_centroid >= r_start & vet_dist_centroid < r_end])
      centroids_df$residue_key <- paste(centroids_df$resid, centroids_df$resno, centroids_df$chain, sep = "_")
      df_selected <- centroids_df[centroids_df$residue_key %in% selected_residues, ]
      
      df_selected_antibody_df <- df_selected[df_selected$chain %in% chain_HL, ]
      df_selected_ligand_df <- df_selected[df_selected$chain %in% chain_AG, ]
      
      local_ab_counts[[step]] <- local_ab_counts[[step]] + table(factor(df_selected_antibody_df$resid, levels = amino_acids))
      local_ag_counts[[step]] <- local_ag_counts[[step]] + table(factor(df_selected_ligand_df$resid, levels = amino_acids))
    }
    
    # --- Return (MODIFICATO: con ok = TRUE) ---
    return(list(
      ok = TRUE,
      matrices_asym = local_matrices_asym, 
      matrices_sym = local_matrices_sym,
      ab_counts = local_ab_counts,
      ag_counts = local_ag_counts
    ))
    
  }, error = function(e) {
    # MODIFICATO: return con 'ok = FALSE' e log dettagliato
    return(list(ok = FALSE, error = paste("Errore PDB:", pdb_path, "Msg:", e$message)))
  })
}


# --- 6. Ciclo principale, aggregazione e salvataggio (MODIFICATO) ---
# Itera su ogni file "group_...txt"
for (list_file_path in list_files_to_process) {
  
  cat("Processando gruppo RADIALE (da lista):", list_file_path, "\n")
  
  # Leggi la lista di PDB per questo gruppo
  pdb_files <- readLines(list_file_path)
  pdb_files <- pdb_files[pdb_files != ""] # Rimuovi righe vuote
  
  if (length(pdb_files) == 0) {
    cat("   ... Lista vuota, salto.\n")
    next
  }
  
  # Esegui l'analisi radiale in parallelo
  results_list <- future_map(pdb_files, process_pdb_radial, .progress = TRUE)
  
  # Separa successi e fallimenti
  valid_results <- keep(results_list, ~ .x$ok)
  failed <- keep(results_list, ~ !.x$ok)
  
  # --- Creazione cartella di output (come Blocco 2) ---
  relative_path <- gsub(paste0("^", list_base_dir, "/"), "", list_file_path)
  output_rel_dir <- gsub("\\.txt$", "", relative_path)
  out_dir <- file.path(radial_matrix_output_dir, output_rel_dir) # Usa la nuova dir di output
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Salva log degli errori
  if (length(failed) > 0) {
    error_messages <- map_chr(failed, "error")
    writeLines(error_messages, file.path(out_dir, "errors.log"))
  }
  
  if (length(valid_results) == 0) {
    cat("   ... Nessun risultato RADIALE valido. Creato file di log.\n")
    next # Salta al prossimo gruppo
  }
  
  # --- INIZIO: Aggregazione specifica per questo gruppo ---
  
  # 1. Inizializza le strutture SOMMA per questo GRUPPO
  group_sum_asym <- setNames(vector("list", length(ring_pairs)), ring_pairs)
  group_sum_sym <- setNames(vector("list", length(ring_pairs)), ring_pairs)
  for (pair in ring_pairs) {
    group_sum_asym[[pair]] <- base_matrix_20x20
    group_sum_sym[[pair]] <- base_matrix_20x20
  }
  
  group_sum_ab_counts <- setNames(vector("list", Nsteps), 1:Nsteps)
  group_sum_ag_counts <- setNames(vector("list", Nsteps), 1:Nsteps)
  for (step in 1:Nsteps) {
    group_sum_ab_counts[[step]] <- base_vector_20
    group_sum_ag_counts[[step]] <- base_vector_20
  }
  
  # 2. Aggrega i risultati (somma i risultati di ogni PDB valido)
  for (res in valid_results) {
    for (pair in ring_pairs) {
      group_sum_asym[[pair]] <- group_sum_asym[[pair]] + res$matrices_asym[[pair]]
      group_sum_sym[[pair]] <- group_sum_sym[[pair]] + res$matrices_sym[[pair]]
    }
    for (step in 1:Nsteps) {
      group_sum_ab_counts[[step]] <- group_sum_ab_counts[[step]] + res$ab_counts[[step]]
      group_sum_ag_counts[[step]] <- group_sum_ag_counts[[step]] + res$ag_counts[[step]]
    }
  }
  
  # --- FINE: Aggregazione ---
  
  # --- Salvataggio (logica dallo snippet, ma dentro il loop) ---
  
  # Salva gli oggetti RDS (liste complete per il gruppo)
  saveRDS(group_sum_asym, file = file.path(out_dir, "contact_matrices_asym_rings.rds"))
  saveRDS(group_sum_sym, file = file.path(out_dir, "contact_matrices_sym_rings.rds"))
  saveRDS(group_sum_ab_counts, file = file.path(out_dir, "residue_counts_interface_ab.rds"))
  saveRDS(group_sum_ag_counts, file = file.path(out_dir, "residue_counts_interface_ligand_df.rds"))
  
  # Salva i CSV individuali per una facile ispezione
  for (pair in ring_pairs) {
    write.csv(group_sum_asym[[pair]], file = file.path(out_dir, paste0("contact_matrices_asym_rings_", pair, ".csv")), row.names = TRUE)
    write.csv(group_sum_sym[[pair]], file = file.path(out_dir, paste0("contact_matrices_sym_rings_", pair, ".csv")), row.names = TRUE)
  }
  for (step in 1:Nsteps) {
    write.csv(group_sum_ab_counts[[step]], file = file.path(out_dir, paste0("residue_counts_interface_ab_step", step, ".csv")), row.names = TRUE)
    write.csv(group_sum_ag_counts[[step]], file = file.path(out_dir, paste0("residue_counts_interface_ligand_df_step", step, ".csv")), row.names = TRUE)
  }
  
  cat("   ... Salvataggio RADIALE completato in:", out_dir, "\n")
}

cat("--- BLOCCO 3: Finito! ---\n")

################################################################################
### QUARTO BLOCCO DI CODICE: CALCOLARE I POTENZIALI STATISTICI RADIALI ###
################################################################################

# --- 1. Definiamo le cartelle di IN e OUT per questo blocco ---

# Input: La cartella di output del Blocco 3, che contiene TUTTI i gruppi
radial_input_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali/matrici_stratificate/"

# Output: Una NUOVA cartella per i potenziali statistici
potential_output_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali/potenziali_statistici_stratificati/"
dir.create(potential_output_dir, showWarnings = FALSE)

# --- 2. Parametri globali e funzioni (dallo snippet) ---
Nsteps <- 3
S <- 0.02
# Lista delle coppie di anelli da analizzare
ring_pairs_list <- list(c(1,1), c(2,2), c(3,3), c(1,2), c(2,3), c(1,3))

# --- Funzione per aggregare CSV SIMMETRICI (dallo snippet) ---
to_sym_df <- function(files) {
  df_list <- list()
  for (i in seq_along(files)) {
    mat <- as.matrix(read.csv(files[i], row.names = 1))
    idx <- which(lower.tri(mat, diag = TRUE), arr.ind = TRUE)
    
    aa1  <- rownames(mat)[idx[, 1]]
    aa2  <- colnames(mat)[idx[, 2]]
    pair <- paste(pmin(aa1, aa2), pmax(aa1, aa2), sep = "-")
    
    pair_name <- sub("\\.csv$", "", basename(files[i]))
    pair_name <- sub("^V_(sym|asym)_", "", pair_name)
    
    df_list[[i]] <- tibble(
      pair = pair,
      !! paste0("V", pair_name) := as.numeric(mat[idx])
    )
  }
  reduce(df_list, inner_join, by = "pair")
}

# --- Funzione per aggregare CSV ASIMMETRICI (dallo snippet) ---
to_asym_df <- function(files) {
  df_list <- list()
  for (i in seq_along(files)) {
    mat <- as.matrix(read.csv(files[i], row.names = 1))
    
    aa1  <- rep(rownames(mat), times = ncol(mat))
    aa2  <- rep(colnames(mat), each  = nrow(mat))
    pair <- paste(aa1, "→", aa2) # Usiamo "→" come nello snippet
    
    pair_name <- sub("\\.csv$", "", basename(files[i]))
    pair_name <- sub("^V_(sym|asym)_", "", pair_name)
    
    df_list[[i]] <- tibble(
      pair = pair,
      !! paste0("V", pair_name) := as.vector(mat)
    )
  }
  reduce(df_list, inner_join, by = "pair")
}

# --- 3. Troviamo tutti i gruppi da processare ---
# Cerchiamo un file "marker" che il Blocco 3 crea in ogni cartella di output
group_marker_files <- list.files(radial_input_dir, 
                                 pattern = "contact_matrices_sym_rings.rds", 
                                 recursive = TRUE, 
                                 full.names = TRUE)

cat("\n--- BLOCCO 4 INIZIATO ---\n")
cat("Trovati", length(group_marker_files), "gruppi da processare per calcolo POTENZIALI.\n")


# --- 4. Ciclo principale: 1 iterazione per GRUPPO ---
for (marker_file in group_marker_files) {
  
  # Definiamo le cartelle IN e OUT per questo specifico gruppo
  group_input_dir <- dirname(marker_file)
  
  # Ricreiamo la stessa struttura di cartelle nell'output dei potenziali
  relative_path <- gsub(paste0("^", radial_input_dir, "/?"), "", group_input_dir)
  group_output_dir <- file.path(potential_output_dir, relative_path)
  dir.create(group_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("Processando potenziali per:", group_input_dir, "\n")
  
  tryCatch({
    
    # --- A. Carica i dati per QUESTO gruppo ---
    contact_matrix_sym <- readRDS(file.path(group_input_dir, "contact_matrices_sym_rings.rds"))
    contact_matrix_asym <- readRDS(file.path(group_input_dir, "contact_matrices_asym_rings.rds"))
    residue_counts_interface_ab <- readRDS(file.path(group_input_dir, "residue_counts_interface_ab.rds"))
    residue_counts_interface_ligando <- readRDS(file.path(group_input_dir, "residue_counts_interface_ligand_df.rds"))
    
    # --- B. Calcola frequenze marginali (Snippet Step 1) ---
    paratope_freq <- list()
    epitope_freq  <- list()
    residue_freq  <- list()
    
    for (step in 1:Nsteps) {
      sum_ab <- sum(residue_counts_interface_ab[[step]])
      sum_lig <- sum(residue_counts_interface_ligando[[step]])
      combined <- residue_counts_interface_ab[[step]] + residue_counts_interface_ligando[[step]]
      sum_comb <- sum(combined)
      
      # Aggiungiamo un controllo per evitare divisioni per zero se un anello è vuoto
      paratope_freq[[step]] <- if (sum_ab > 0) residue_counts_interface_ab[[step]] / sum_ab else residue_counts_interface_ab[[step]]
      epitope_freq[[step]]  <- if (sum_lig > 0) residue_counts_interface_ligando[[step]] / sum_lig else residue_counts_interface_ligando[[step]]
      residue_freq[[step]] <- if (sum_comb > 0) combined / sum_comb else combined
    }
    
    # --- C. Calcola potenziali per coppia di anelli (Snippet Step 2) ---
    for (idx in seq_along(ring_pairs_list)) {
      r1 <- ring_pairs_list[[idx]][1]
      r2 <- ring_pairs_list[[idx]][2]
      pair_name <- paste0(r1, "-", r2)
      
      P_par  <- paratope_freq[[r1]]
      P_epi  <- epitope_freq[[r2]]
      P_sym1 <- residue_freq[[r1]]
      P_sym2 <- residue_freq[[r2]]
      
      cm_asym <- contact_matrix_asym[[pair_name]]
      cm_sym  <- contact_matrix_sym[[pair_name]]
      
      # Frequenze osservate (con controllo per divisione per zero)
      sum_cm_asym <- sum(cm_asym)
      P_obs_asym <- if (sum_cm_asym > 0) cm_asym / sum_cm_asym else cm_asym
      
      sum_cm_sym <- sum(cm_sym[lower.tri(cm_sym, diag = TRUE)])
      P_obs_sym  <- if (sum_cm_sym > 0) cm_sym / sum_cm_sym else cm_sym
      
      # Frequenze attese
      outer_asym <- outer(P_par, P_epi, "*")
      outer_sym  <- outer(P_sym1, P_sym2, "*")
      
      # Calcolo Ratio (con gestione di Inf)
      ratio_asym <- P_obs_asym / outer_asym
      ratio_sym  <- P_obs_sym  / outer_sym
      
      # Calcolo Potenziali V (formula dallo snippet)
      V_asym <- (log(1 + cm_asym * S) - log(1 + cm_asym * S * ratio_asym)) * 2.479
      V_sym  <- (log(1 + cm_sym  * S) - log(1 + cm_sym  * S * ratio_sym)) * 2.479
      
      # Pulizia (dallo snippet)
      V_asym[!is.finite(V_asym)] <- 0
      V_sym[!is.finite(V_sym)]   <- 0
      
      # --- D. Salva i file CSV intermedi NELLA cartella del gruppo ---
      write.csv(V_asym, file.path(group_output_dir, paste0("V_asym_", r1, "_", r2, ".csv")))
      write.csv(V_sym,  file.path(group_output_dir, paste0("V_sym_", r1, "_", r2, ".csv")))
    }
    
    # --- E. Trova i CSV appena creati ---
    sym_files  <- list.files(group_output_dir, pattern = "^V_sym_.*\\.csv$", full.names = TRUE)
    asym_files <- list.files(group_output_dir, pattern = "^V_asym_.*\\.csv$", full.names = TRUE)
    
    # Ordina i file per coerenza
    sym_files  <- sym_files[order(sym_files)]
    asym_files <- asym_files[order(asym_files)]
    
    # --- F. Aggrega in DataFrame Wide e salva ---
    if (length(sym_files) == 6) { # Assicurati che ci siano tutti e 6
      df_sym_wide  <- to_sym_df(sym_files)
      df_sym_wide <- df_sym_wide %>%
        separate(pair, into = c("aa1", "aa2"), sep = "-") %>%
        relocate(aa1, aa2)
      
      write.csv(df_sym_wide, file.path(group_output_dir, "potentials_sym_wide.csv"), row.names = FALSE)
      saveRDS(df_sym_wide, file.path(group_output_dir, "potentials_sym_wide.rds"))
    }
    
    if (length(asym_files) == 6) { # Assicurati che ci siano tutti e 6
      df_asym_wide <- to_asym_df(asym_files)
      # Completo lo snippet che era interrotto
      df_asym_wide <- df_asym_wide %>%
        separate(pair, into = c("aa1", "aa2"), sep = " → ") %>%
        relocate(aa1, aa2)
      
      write.csv(df_asym_wide, file.path(group_output_dir, "potentials_asym_wide.csv"), row.names = FALSE)
      saveRDS(df_asym_wide, file.path(group_output_dir, "potentials_asym_wide.rds"))
    }
    
    cat("   ... Potenziali salvati in:", group_output_dir, "\n")
    
  }, error = function(e) {
    # Se un gruppo fallisce, logga l'errore e continua con il prossimo
    error_msg <- paste("!!! ERRORE processando gruppo:", group_input_dir, "Msg:", e$message)
    cat(error_msg, "\n")
    writeLines(error_msg, file.path(group_output_dir, "error.log"))
  })
}

cat("--- BLOCCO 4: Finito! ---\n")

### Alla fine di tutto ciò ho 3740 potenziali statistici in tutto (che sono chiaramente troppi)
### Non posso utilizzare tutti e 3740 i potenziali per assegnare i punteggi alle 22000 pose di docking poiché ciò richiederebbe circa 2 ore per potenziale.
### Decido quindi di seguire un altro approccio. Metto insieme tutti i potenziali di tutte le seed e per ogni group (group_2, group_5 etc...) scelgo randomicamente 10 potenziali simmetrici e 10 potenziali asimmetrici

################################################################################
### QUINTO BLOCCO DI CODICE: CAMPIONAMENTO DEI POTENZIALI ###
################################################################################

# --- 1. Definiamo le cartelle di IN e OUT per questo blocco ---

# Input: La cartella di output del Blocco 4, che contiene TUTTI i potenziali
potentials_input_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali/potenziali_statistici_stratificati/"

# Output: Una NUOVA cartella per le LISTE di file .txt campionati
selected_potentials_list_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali/potenziali_selezionati_liste/"
dir.create(selected_potentials_list_dir, showWarnings = FALSE)

# --- 2. Impostiamo un seed per la riproducibilità ---
# Questo è FONDAMENTALE per poter rieseguire e ottenere lo STESSO campionamento
set.seed(12345) 

cat("\n--- BLOCCO 5 INIZIATO ---\n")
cat("Campionamento casuale dei potenziali...\n")

# --- 3. Troviamo tutti i file .rds dei potenziali ---
all_sym_files <- list.files(potentials_input_dir, 
                            pattern = "potentials_sym_wide.rds", 
                            recursive = TRUE, 
                            full.names = TRUE)

all_asym_files <- list.files(potentials_input_dir, 
                             pattern = "potentials_asym_wide.rds", 
                             recursive = TRUE, 
                             full.names = TRUE)

cat("Trovati", length(all_sym_files), "potenziali simmetrici e", length(all_asym_files), "asimmetrici.\n")

# --- 4. Creiamo un database di percorsi ed estraiamo il gruppo 'split_n' ---
sym_potentials_db <- tibble(path = all_sym_files) %>%
  mutate(split_group = stringr::str_extract(path, "split_\\d+"))

asym_potentials_db <- tibble(path = all_asym_files) %>%
  mutate(split_group = stringr::str_extract(path, "split_\\d+"))

# --- 5. Eseguiamo il campionamento (10 per ogni 'split_group') ---

# SIMMETRICI
selected_sym_lists <- sym_potentials_db %>%
  group_by(split_group) %>%
  sample_n(2, replace = FALSE) %>% # 'replace = FALSE' assicura che siano 10 unici
  summarise(paths = list(path), .groups = 'drop') # Raggruppa i 10 path in una lista

# ASIMMETRICI
selected_asym_lists <- asym_potentials_db %>%
  group_by(split_group) %>%
  sample_n(2, replace = FALSE) %>%
  summarise(paths = list(path), .groups = 'drop')

# --- 6. Salviamo le liste di percorsi in file .txt ---

# Salvataggio liste simmetriche
for (i in 1:nrow(selected_sym_lists)) {
  split_name <- selected_sym_lists$split_group[i]
  paths_to_save <- unlist(selected_sym_lists$paths[i])
  output_file <- file.path(selected_potentials_list_dir, paste0("selected_sym_", split_name, ".txt"))
  
  writeLines(paths_to_save, con = output_file)
}
cat("Salvate", nrow(selected_sym_lists), "liste di potenziali SIMMETRICI campionati in:\n", selected_potentials_list_dir, "\n")

# Salvataggio liste asimmetriche
for (i in 1:nrow(selected_asym_lists)) {
  split_name <- selected_asym_lists$split_group[i]
  paths_to_save <- unlist(selected_asym_lists$paths[i])
  output_file <- file.path(selected_potentials_list_dir, paste0("selected_asym_", split_name, ".txt"))
  
  writeLines(paths_to_save, con = output_file)
}
cat("Salvate", nrow(selected_asym_lists), "liste di potenziali ASIMMETRICI campionati in:\n", selected_potentials_list_dir, "\n")


cat("--- BLOCCO 5: Finito! ---\n")

# PROCEDO A SCRIVERE IL NUOVO SCRIPT CHE LANCERO' CON CAFFEINATE NOHUP