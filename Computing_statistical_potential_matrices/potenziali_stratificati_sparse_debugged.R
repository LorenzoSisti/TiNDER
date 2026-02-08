# THIS SCRIPT COMPUTES AMINO-ACID STATISTICAL POTENTIALS WITHIN CONCENTRIC RINGS
# AROUND THE ANTIBODY–ANTIGEN INTERFACE. IT LOADS PRECOMPUTED CONTACT MATRICES
# AND RESIDUE COUNTS, NORMALIZES FREQUENCIES, DERIVES ASYMMETRIC/SYMMETRIC
# POTENTIALS VIA LOG-RATIOS, SAVES CSV OUTPUTS, AND PLOTS HEATMAPS PER RING.

### Libraries ###
library(pheatmap)
library(patchwork)
library(ggplotify)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)


### Setting paths, directories and global variables ###
pdb_dir <- "/Users/lorenzosisti/Downloads/database_settembre_renamed/"
results_dir <- "/Users/lorenzosisti/Downloads/matrici_stratificate_sparse/" 
dir.create(results_dir, showWarnings = FALSE)


DistCutoff <- 8.5
Nsteps <- 3
S <- 0.02

amino_acids <- c("ARG","LYS","ASN","ASP","GLN","GLU","HIS","PRO","TYR","TRP",
                 "SER","THR","GLY","ALA","MET","CYS","PHE","LEU","VAL","ILE")

### Load saved data ###
contact_matrix_sym <- readRDS(file.path(results_dir, "contact_matrices_sym_rings.rds"))
contact_matrix_asym <- readRDS(file.path(results_dir, "contact_matrices_asym_rings.rds"))
residue_counts_interface_ab <- readRDS(file.path(results_dir, "residue_counts_interface_ab.rds"))
residue_counts_interface_ligando <- readRDS(file.path(results_dir, "residue_counts_interface_ligand_df.rds"))


### === STEP 1: compute marginal frequencies per ring === ###
paratope_freq <- list()
epitope_freq  <- list()
residue_freq  <- list()

for (step in 1:Nsteps) {
  
  paratope_freq[[step]] <- residue_counts_interface_ab[[step]] /
    sum(residue_counts_interface_ab[[step]])
  
  epitope_freq[[step]]  <- residue_counts_interface_ligando[[step]] /
    sum(residue_counts_interface_ligando[[step]])
  
  combined <- residue_counts_interface_ab[[step]] + residue_counts_interface_ligando[[step]]
  residue_freq[[step]] <- combined / sum(combined)
}

### === STEP 2: compute potentials for each pair of rings === ###
ring_pairs <- list(c(1,1), c(2,2), c(3,3), c(1,2), c(2,3), c(1,3))

for (idx in seq_along(ring_pairs)) {
  
  r1 <- ring_pairs[[idx]][1]
  r2 <- ring_pairs[[idx]][2]
  
  pair_name <- paste0(r1, "-", r2)
  
  P_par  <- paratope_freq[[r1]]
  P_epi  <- epitope_freq[[r2]]
  P_sym1 <- residue_freq[[r1]]
  P_sym2 <- residue_freq[[r2]]
  
  cm_asym <- contact_matrix_asym[[pair_name]]
  cm_sym  <- contact_matrix_sym[[pair_name]]
  
  P_obs_asym <- cm_asym / sum(cm_asym)
  P_obs_sym  <- cm_sym  / sum(cm_sym[lower.tri(cm_sym, diag = TRUE)])
  
  ### ✅ NUOVA FORMULA: Queste sono matrici relative
  ratio_asym <- P_obs_asym / outer(P_par, P_epi, "*")
  ratio_sym  <- P_obs_sym  / outer(P_sym1, P_sym2, "*")
  
  V_asym <- (log(1 + cm_asym * S) - log(1 + cm_asym * S * ratio_asym)) * 2.479
  V_sym  <- (log(1 + cm_sym  * S) - log(1 + cm_sym  * S * ratio_sym)) * 2.479
  
  ### pulizia infinities / NaN
  V_asym[!is.finite(V_asym)] <- 0
  V_sym[!is.finite(V_sym)]   <- 0
  
  write.csv(V_asym, file.path(results_dir, paste0("V_asym_", r1, "_", r2, ".csv")))
  write.csv(V_sym,  file.path(results_dir, paste0("V_sym_", r1, "_", r2, ".csv")))
}


# 1. Lista dei nomi file prodotti al passo precedente
sym_files  <- list.files(results_dir, pattern = "^V_sym_.*\\.csv$", full.names = TRUE)
asym_files <- list.files(results_dir, pattern = "^V_asym_.*\\.csv$", full.names = TRUE)

# ✅ Ordina per numerazione (1-1, 1-2, 2-2, ...)
sym_files  <- sym_files[order(sym_files)]
asym_files <- asym_files[order(asym_files)]


### ---- FUNZIONE PER SIMMETRICI (210 righe) ---- ###
to_sym_df <- function(files) {
  
  df_list <- list()
  
  for (i in seq_along(files)) {
    
    mat <- as.matrix(read.csv(files[i], row.names = 1))
    idx <- which(lower.tri(mat, diag = TRUE), arr.ind = TRUE)
    
    aa1  <- rownames(mat)[idx[, 1]]
    aa2  <- colnames(mat)[idx[, 2]]
    pair <- paste(pmin(aa1, aa2), pmax(aa1, aa2), sep = "-")
    
    # ✅ nome colonna: 1_1 , 1_2 , 1_3 , 2_2 , ...
    pair_name <- sub("\\.csv$", "", basename(files[i]))
    pair_name <- sub("^V_(sym|asym)_", "", pair_name)
    
    df_list[[i]] <- tibble(
      pair = pair,
      !! paste0("V", pair_name) := as.numeric(mat[idx])
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
    
    # ✅ stesso fix per l'asimmetrico
    pair_name <- sub("\\.csv$", "", basename(files[i]))
    pair_name <- sub("^V_(sym|asym)_", "", pair_name)
    
    df_list[[i]] <- tibble(
      pair = pair,
      !! paste0("V", pair_name) := as.vector(mat)
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
  relocate(aa1, aa2)   # porta le nuove colonne davanti

# ASIMMETRICO (pair ha formato "ARG → GLY")
df_asym_wide <- df_asym_wide %>%
  separate(pair, into = c("aa1", "aa2"), sep = " → ") %>%
  relocate(aa1, aa2)

# Aggiunge il suffisso alle colonne aa1 e aa2 SOLO per df_asym_wide
df_asym_wide <- df_asym_wide %>%
  mutate(
    aa1 = paste0(aa1, "_Ab"),
    aa2 = paste0(aa2, "_Ag")
  )



### ✅ Salvataggio finale
write.csv(df_sym_wide,  file.path(results_dir, "V_sym_wide.csv"),  row.names = FALSE)
write.csv(df_asym_wide, file.path(results_dir, "V_asym_wide.csv"), row.names = FALSE)

message("\n✅ Generati:\n- V_sym_wide.csv (210 righe × 6 colonne corsi: V1-1, V2-2, V3-3, V1-2, V2-3, V1-3)\n- V_asym_wide.csv (400 righe × 6 colonne corsi: V1-1, V2-2, V3-3, V1-2, V2-3, V1-3)\n")
