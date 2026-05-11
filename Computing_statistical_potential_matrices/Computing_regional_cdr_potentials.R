# Caricamento delle librerie necessarie
library(bio3d)  
library(reshape2)  
library(ftrCOOL)  
library(AATtools)  
library(ggplot2)
library(tidyr)
library(purrr)
library(dplyr)    # Necessario per %>% e group_by()

# Cartella contenente i file PDB
pdb_folder <- "/Users/lorenzosisti/Downloads/database_settembre_renamed/"
results_dir <- "/Users/lorenzosisti/Downloads/potenziali_statistici_cdr/"
dir.create(results_dir, showWarnings = FALSE)

pdb_files <- list.files(path = pdb_folder, pattern = "\\.pdb$", full.names = TRUE)

# Cutoff per il contatto (in Å)  
DistCutoff <- 8.5
S <- 0.02

# Amino acids ordered according to the Kyte-Doolittle hydrophobicity scale
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE") 

set.seed(1234)

# Inizializzazione delle strutture dati cumulative
cdr_index <- c("cdr_h1", "cdr_h2", "cdr_h3", "cdr_l1", "cdr_l2", "cdr_l3")

# Definizione strutturata delle CDR per facilitare l'iterazione nel ciclo for
cdr_definitions <- list(
  cdr_h1 = list(chain = "H", resno = 26:32),
  cdr_h2 = list(chain = "H", resno = 52:56),
  cdr_h3 = list(chain = "H", resno = 95:102),
  cdr_l1 = list(chain = "L", resno = 24:34),
  cdr_l2 = list(chain = "L", resno = 50:56),
  cdr_l3 = list(chain = "L", resno = 89:97)
)

# Inizializzamo le strutture dati CUMULATIVE per le 6 matrici
contact_matrices_asym_rings_sum <- vector("list", length(cdr_index))
names(contact_matrices_asym_rings_sum) <- cdr_index

contact_matrices_sym_rings_sum  <- vector("list", length(cdr_index))
names(contact_matrices_sym_rings_sum) <- cdr_index

residue_counts_interface_ab_sum <- vector("list", length(cdr_index)) 
names(residue_counts_interface_ab_sum) <- cdr_index

# Residui ligando (antigene) cumulativi (uno solo globale per l'intera interfaccia)
residue_counts_interface_ligand_sum <- setNames(rep(0, length(amino_acids)), amino_acids)

# Definiamo una matrice 20x20 quadrata da inizializzare a zero
base_matrix <- matrix(0, nrow = 20, ncol = 20, dimnames = list(amino_acids, amino_acids))

for (index in cdr_index) {
  contact_matrices_asym_rings_sum[[index]] <- base_matrix
  contact_matrices_sym_rings_sum[[index]]  <- base_matrix
  residue_counts_interface_ab_sum[[index]] <- setNames(rep(0, length(amino_acids)), amino_acids)
}

#pdb_file <- "/Users/lorenzosisti/Downloads/database_settembre_renamed//1afv.pdb"

# Inizio Loop su tutti i file PDB  
for (pdb_file in pdb_files) { 
  
  tryCatch({
    
    # Lettura della struttura PDB  
    pdb_aus <- read.pdb(pdb_file) 
    chain_HL <- c("H", "L")
    chain_AG <- "A" # Assicurati che l'antigene sia sempre "A" oppure usa un'esclusione come nel tuo primo script
    
    # Coarse-grained representation
    df_coord <- pdb_aus$atom[pdb_aus$atom$type != "HETATM" & (
      (pdb_aus$atom$resid != "GLY" & !(pdb_aus$atom$elety %in% c("N", "CA", "C", "O"))) |
        (pdb_aus$atom$resid == "GLY" & pdb_aus$atom$elety == "CA")),]
    
    if (nrow(df_coord) == 0) {
      cat("File con df_coord vuoto:", pdb_file, "\n", file = "errors.log", append = TRUE)
      next # Salta al prossimo PDB
    }
    
    df_coord_corrected <- df_coord
    
    # 1. Gestione dei valori NA nella colonna insert (li trasformiamo in stringhe vuote)
    df_coord_corrected$insert[is.na(df_coord_corrected$insert)] <- ""
    
    # 2. Compute centroid coordinates (Aggiunto 'insert' nel group_by)
    centroidi_df <- as.data.frame(
      df_coord_corrected %>%
        group_by(chain, resno, insert, resid) %>%
        summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
    )
    
    # 3. Creazione di nomi univoci che includono la lettera di inserzione (es. PHE_100A_H)
    res_names <- paste(centroidi_df$resid, 
                       paste0(centroidi_df$resno, centroidi_df$insert), 
                       centroidi_df$chain, sep = "_")
    
    rownames(centroidi_df) <- res_names
    
    # 4. Salvataggio dataframe con le coordinate per la matrice di distanza
    df_coord_resid_xyz <- centroidi_df[, c("resid", "resno", "x", "y", "z")] 
    rownames(df_coord_resid_xyz) <- res_names
    
    # Compute Distance Matrix
    DistMat <- as.matrix(dist(df_coord_resid_xyz[, c("x", "y", "z")])) 
    DistMat_bin <- ifelse(DistMat <= DistCutoff, 1, 0) 
    
    # Define logical conditions
    condA <- !(centroidi_df$chain %in% c("H", "L")) # Antigene (o usa centroidi_df$chain == chain_AG)
    condHL <- centroidi_df$chain %in% c("H", "L")   # Anticorpo
    
    # Subset distance matrix (Righe = Antigene, Colonne = Anticorpo)
    Inter_DistMat_Bin <- DistMat_bin[condA, condHL, drop = FALSE]
    
    # Vettore per tracciare tutti i residui dell'antigene che toccano una qualsiasi CDR
    ag_interacting_with_any_cdr <- c()
    
    # =========================================================
    # ITERAZIONE SULLE 6 CDR
    # =========================================================
    for (cdr_name in cdr_index) {
      
      # Estrazione info specifiche per questa CDR
      cdr_chain <- cdr_definitions[[cdr_name]]$chain
      cdr_resnos <- cdr_definitions[[cdr_name]]$resno
      
      # Identifica i nomi dei residui di QUESTA CDR
      cond_this_cdr <- centroidi_df$chain == cdr_chain & centroidi_df$resno %in% cdr_resnos
      this_cdr_residue_names <- rownames(centroidi_df)[cond_this_cdr]
      
      # Mantiene solo le colonne della matrice corrispondenti a questa CDR
      valid_cols <- intersect(colnames(Inter_DistMat_Bin), this_cdr_residue_names)
      
      if (length(valid_cols) == 0) {
        next # Salta se non ci sono residui di questa CDR nella struttura
      }
      
      Inter_DistMat_Bin_CDR <- Inter_DistMat_Bin[, valid_cols, drop = FALSE]
      
      # Calcola i contatti specifici per questa CDR
      BS_A_cdr <- apply(Inter_DistMat_Bin_CDR, 1, sum)
      BS_A_cdr <- BS_A_cdr[BS_A_cdr != 0]
      
      BS_HL_cdr <- apply(Inter_DistMat_Bin_CDR, 2, sum)
      BS_HL_cdr <- BS_HL_cdr[BS_HL_cdr != 0]
      
      BS_A_names_cdr <- names(BS_A_cdr)
      BS_HL_names_cdr <- names(BS_HL_cdr)
      
      # Aggiunge i residui Ag all'insieme globale di questo PDB
      ag_interacting_with_any_cdr <- unique(c(ag_interacting_with_any_cdr, BS_A_names_cdr))
      
      # Aggiorna le matrici cumulative per questa CDR
      for (res_ab in BS_HL_names_cdr) {
        idx_ab <- which(colnames(Inter_DistMat_Bin_CDR) == res_ab)
        contacts <- which(Inter_DistMat_Bin_CDR[, idx_ab] == 1)
        
        for (idx_ag in contacts) {
          res_ag <- rownames(Inter_DistMat_Bin_CDR)[idx_ag]
          aa_ab <- strsplit(res_ab, "_")[[1]][1]
          aa_ag <- strsplit(res_ag, "_")[[1]][1]
          
          if (aa_ab %in% amino_acids && aa_ag %in% amino_acids) {
            
            # Matrice Asimmetrica
            contact_matrices_asym_rings_sum[[cdr_name]][aa_ab, aa_ag] <- contact_matrices_asym_rings_sum[[cdr_name]][aa_ab, aa_ag] + 1
            
            # Matrice Simmetrica
            if (aa_ab == aa_ag) {
              contact_matrices_sym_rings_sum[[cdr_name]][aa_ab, aa_ag] <- contact_matrices_sym_rings_sum[[cdr_name]][aa_ab, aa_ag] + 1
            } else {
              contact_matrices_sym_rings_sum[[cdr_name]][aa_ab, aa_ag] <- contact_matrices_sym_rings_sum[[cdr_name]][aa_ab, aa_ag] + 1
              contact_matrices_sym_rings_sum[[cdr_name]][aa_ag, aa_ab] <- contact_matrices_sym_rings_sum[[cdr_name]][aa_ag, aa_ab] + 1
            }
          }
        }
      }
      
      # Conteggio dei residui per questa CDR
      df_centroidi_cdr_bs <- centroidi_df[rownames(centroidi_df) %in% BS_HL_names_cdr, ]
      if (nrow(df_centroidi_cdr_bs) > 0) {
        residue_counts_interface_ab_sum[[cdr_name]] <- residue_counts_interface_ab_sum[[cdr_name]] + table(factor(df_centroidi_cdr_bs$resid, levels = amino_acids))
      }
      
    } # Fine iterazione CDR
    
    # Aggiornamento residui antigene (una sola volta per PDB, usando la lista unica dei contatti)
    if (length(ag_interacting_with_any_cdr) > 0) {
      df_centroidi_ag_bs <- centroidi_df[rownames(centroidi_df) %in% ag_interacting_with_any_cdr, ]
      residue_counts_interface_ligand_sum <- residue_counts_interface_ligand_sum + table(factor(df_centroidi_ag_bs$resid, levels = amino_acids))
    }
    
  }, error = function(e) {
    # Logging degli errori senza interrompere l'intero loop
    cat(paste("Errore elaborazione file:", pdb_file, "-", e$message, "\n"), file = "errors.log", append = TRUE)
  })
  
}

# =========================================================
# SALVATAGGIO DEI RISULTATI IN RDS E CSV
# =========================================================

cat("Inizio salvataggio dei risultati in:", results_dir, "\n")

# 1. SALVATAGGIO IN FORMATO RDS (Mantiene la struttura ad albero/lista di R)
saveRDS(contact_matrices_asym_rings_sum, file = file.path(results_dir, "contact_matrices_asym_all_cdrs.rds"))
saveRDS(contact_matrices_sym_rings_sum, file = file.path(results_dir, "contact_matrices_sym_all_cdrs.rds"))
saveRDS(residue_counts_interface_ab_sum, file = file.path(results_dir, "residue_counts_ab_all_cdrs.rds"))
saveRDS(residue_counts_interface_ligand_sum, file = file.path(results_dir, "residue_counts_ag_global.rds"))

# 2. SALVATAGGIO IN FORMATO CSV (Esporta un file per ogni CDR per leggibilità)

# Ciclo sulle 6 CDR per spacchettare le liste
for (cdr_name in cdr_index) {
  
  # Salva Matrice Asimmetrica
  file_asym <- file.path(results_dir, paste0("Contact_Matrix_Asym_", cdr_name, ".csv"))
  write.csv(contact_matrices_asym_rings_sum[[cdr_name]], file = file_asym, row.names = TRUE)
  
  # Salva Matrice Simmetrica
  file_sym <- file.path(results_dir, paste0("Contact_Matrix_Sym_", cdr_name, ".csv"))
  write.csv(contact_matrices_sym_rings_sum[[cdr_name]], file = file_sym, row.names = TRUE)
  
  # Salva Conteggi Residui Anticorpo (Convertiamo la table in data.frame per un CSV più pulito)
  file_counts_ab <- file.path(results_dir, paste0("Residue_Counts_Ab_", cdr_name, ".csv"))
  df_counts_ab <- as.data.frame(residue_counts_interface_ab_sum[[cdr_name]])
  colnames(df_counts_ab) <- c("AminoAcid", "Count")
  write.csv(df_counts_ab, file = file_counts_ab, row.names = FALSE)
}

# Salva Conteggi Residui Antigene (Questo è globale, fuori dal loop delle CDR)
file_counts_ag <- file.path(results_dir, "Residue_Counts_Ag_Global.csv")
df_counts_ag <- as.data.frame(residue_counts_interface_ligand_sum)
colnames(df_counts_ag) <- c("AminoAcid", "Count")
write.csv(df_counts_ag, file = file_counts_ag, row.names = FALSE)

cat("Salvataggio completato con successo!\n")

### === STEP 1: Calcolo frequenze marginali per CDR e Antigene === ###
cdr_freq <- list()
sym_freq <- list() # Per le matrici simmetriche

# La frequenza dell'antigene è globale per l'interfaccia
antigen_freq <- residue_counts_interface_ligand_sum / sum(residue_counts_interface_ligand_sum)

for (cdr_name in cdr_index) {
  
  # Frequenze per l'anticorpo (specifiche per ogni CDR)
  ab_counts <- residue_counts_interface_ab_sum[[cdr_name]]
  cdr_freq[[cdr_name]] <- ab_counts / sum(ab_counts)
  
  # Frequenze combinate per il potenziale simmetrico (CDR + Antigene)
  total_counts <- ab_counts + residue_counts_interface_ligand_sum
  sym_freq[[cdr_name]] <- total_counts / sum(total_counts)
}

### === STEP 2: Calcolo Potenziali Statistici (6 coppie: CDR vs Ag) === ###

for (cdr_name in cdr_index) {
  
  # Frequenze marginali per questa interazione
  P_par  <- cdr_freq[[cdr_name]]
  P_epi  <- antigen_freq
  P_sym  <- sym_freq[[cdr_name]] 
  
  cm_asym <- contact_matrices_asym_rings_sum[[cdr_name]]
  cm_sym  <- contact_matrices_sym_rings_sum[[cdr_name]]
  
  # Probabilità osservate
  P_obs_asym <- cm_asym / sum(cm_asym)
  P_obs_sym  <- cm_sym  / sum(cm_sym)
  
  ### Matrici relative (Expected Probabilities)
  ratio_asym <- P_obs_asym / outer(P_par, P_epi, "*")
  ratio_sym  <- P_obs_sym  / outer(P_sym, P_sym, "*") 
  
  ### Formula Log-Odds
  V_asym <- (log(1 + cm_asym * S) - log(1 + cm_asym * S * ratio_asym)) * 2.479
  V_sym  <- (log(1 + cm_sym  * S) - log(1 + cm_sym  * S * ratio_sym)) * 2.479
  
  ### Pulizia infinities / NaN (se cm è 0 o frequenze mancano)
  V_asym[!is.finite(V_asym)] <- 0
  V_sym[!is.finite(V_sym)]   <- 0
  
  # Salvataggio
  write.csv(V_asym, file.path(results_dir, paste0("V_asym_", cdr_name, ".csv")))
  write.csv(V_sym,  file.path(results_dir, paste0("V_sym_", cdr_name, ".csv")))
}

# =========================================================
# STEP 3: POST-PROCESSING & COSTRUZIONE DATAFRAME WIDE
# =========================================================

# Recupero file generati
sym_files  <- list.files(results_dir, pattern = "^V_sym_cdr_.*\\.csv$", full.names = TRUE)
asym_files <- list.files(results_dir, pattern = "^V_asym_cdr_.*\\.csv$", full.names = TRUE)

# Riordino forzato dei file per rispettare l'ordine originale di cdr_index
order_files <- function(files) {
  names_in_files <- sub("\\.csv$", "", sub("^V_(sym|asym)_", "", basename(files)))
  files[match(cdr_index, names_in_files)]
}
sym_files  <- order_files(sym_files)
asym_files <- order_files(asym_files)

### ---- FUNZIONE PER SIMMETRICI (210 righe) ---- ###
to_sym_df <- function(files) {
  
  df_list <- list()
  
  for (i in seq_along(files)) {
    mat <- as.matrix(read.csv(files[i], row.names = 1))
    idx <- which(lower.tri(mat, diag = TRUE), arr.ind = TRUE)
    
    aa1  <- rownames(mat)[idx[, 1]]
    aa2  <- colnames(mat)[idx[, 2]]
    pair <- paste(pmin(aa1, aa2), pmax(aa1, aa2), sep = "-")
    
    # Nome colonna derivato dal file (es. cdr_h1)
    pair_name <- sub("\\.csv$", "", basename(files[i]))
    pair_name <- sub("^V_sym_", "", pair_name) 
    
    df_list[[i]] <- tibble(
      pair = pair,
      !! paste0("V_", pair_name) := as.numeric(mat[idx])
    )
  }
  
  reduce(df_list, inner_join, by = "pair")
}

### ---- FUNZIONE PER ASIMMETRICI (400 righe) ---- ###
to_asym_df <- function(files) {
  
  df_list <- list()
  
  for (i in seq_along(files)) {
    mat <- as.matrix(read.csv(files[i], row.names = 1))
    
    aa1  <- rep(rownames(mat), times = ncol(mat))
    aa2  <- rep(colnames(mat), each  = nrow(mat))
    pair <- paste(aa1, "→", aa2)
    
    pair_name <- sub("\\.csv$", "", basename(files[i]))
    pair_name <- sub("^V_asym_", "", pair_name)
    
    df_list[[i]] <- tibble(
      pair = pair,
      !! paste0("V_", pair_name) := as.vector(mat)
    )
  }
  
  reduce(df_list, inner_join, by = "pair")
}

### ✅ Costruzione DataFrame WIDE
df_sym_wide  <- to_sym_df(sym_files)
df_asym_wide <- to_asym_df(asym_files)

### ---- POST-PROCESSING: estrai aa1 e aa2 dal campo `pair` ---- ###

# SIMMETRICO (pair ha formato "ARG-GLY")
df_sym_wide <- df_sym_wide %>%
  separate(pair, into = c("aa1", "aa2"), sep = "-") %>%
  relocate(aa1, aa2)

# ASIMMETRICO (pair ha formato "ARG → GLY")
df_asym_wide <- df_asym_wide %>%
  separate(pair, into = c("aa1", "aa2"), sep = " → ") %>%
  relocate(aa1, aa2) %>%
  mutate(
    aa1 = paste0(aa1, "_Ab"),
    aa2 = paste0(aa2, "_Ag")
  )

### ✅ Salvataggio finale
write.csv(df_sym_wide,  file.path(results_dir, "V_sym_wide.csv"),  row.names = FALSE)
write.csv(df_asym_wide, file.path(results_dir, "V_asym_wide.csv"), row.names = FALSE)

message("\n✅ Generati:\n- V_sym_wide.csv (210 righe × 8 colonne: aa1, aa2, V_cdr_h1, ..., V_cdr_l3)\n- V_asym_wide.csv (400 righe × 8 colonne: aa1, aa2, V_cdr_h1, ..., V_cdr_l3)\n")


