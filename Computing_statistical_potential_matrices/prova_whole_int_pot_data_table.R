#########################################################################################
# THIS SCRIPT COMPUTES AMINO ACID CONTACT MATRICES AND WHOLE-INTERFACE STATISTICAL POTENTIALS FROM A DIRECTORY OF PDB STRUCTURES. 
# IT ANALYZES ANTIBODY–ANTIGEN INTERACTIONS BY:
# 1) CALCULATING SIDE-CHAIN CENTROIDS AND CONTACT MATRICES
# 2) IDENTIFYING INTERFACE RESIDUES BETWEEN ANTIBODY AND ANTIGEN
# 3) SUMMARIZING CONTACT FREQUENCIES ACROSS ALL STRUCTURES
# 4) COMPUTING ASYMMETRIC AND SYMMETRIC WHOLE-INTERFACE STATISTICAL POTENTIALS
# 5) VISUALIZING THE RESULTS AS HEATMAPS.
#########################################################################################

# DA IMPLEMENTARE FUNZIONE PER POTENZIALI SIMMETRICI
# DA IMPLEMENTARE CALCOLO pmf SECONDO SIPPL

### Required libraries
pacman::p_load(bio3d,dplyr,future,furrr,purrr,progressr,pheatmap,patchwork,ggplotify,reshape2,tidyr,data.table)

# Define the path to a custom function files
#source("/path/to/your/custom/functions/files/functions.R")
source("/Users/lorenzosisti/Documents/Script_ottimizzati_funzioni/functions.R")

### Set up parallelization to speed up computation
plan(multisession, workers = parallel::detectCores() - 1)
handlers(global = TRUE)
handlers("rstudio")

### Define directories and global parameters
pdb_dir <- "/Users/lorenzosisti/Downloads/database_settembre_renamed/"
results_dir <- "/Users/lorenzosisti/Downloads/potenziali_statistici_whole_26_06_data_table_sippl/"
dir.create(results_dir, showWarnings = FALSE)

# Distance cutoff (Å) to define contact between side-chains centroids
DistCutoff <- 8.5  
S <- 0.02

# Amino acids ordered according to the Kyte-Doolittle hydrophobicity scale
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE") 

set.seed(1234)

all_pdbs <- list.files(pdb_dir, pattern = "*.pdb", recursive = TRUE, full.names = TRUE)

# pdb_path <- "/Users/lorenzosisti/Downloads/database_settembre_renamed//9rm2.pdb" 

### Main processing function
gen_df_contacts <- function(pdb_path) {
  aa <- aa.table$aa3[1:20]
  file_name <- basename(pdb_path)
  tryCatch({
    pdb_aus <- read.pdb(pdb_path)
    renumbered_df <- renumber_ab_chains(pdb_aus, pdb_path, log_file = "errors.log")
    if (!renumbered_df$ok) {
      return(list(ok = FALSE, filename = file_name, path = pdb_path, error = renumbered_df$error))
    }
    dt_coord <- as.data.table(renumbered_df$df_coord_renumbered)
    dt_centroids <- dt_coord[, .(x = mean(x), y = mean(y), z = mean(z)), by = .(chain, resno, resid)]
    dt_ab <- dt_centroids[chain %in% c("H", "L")]
    dt_ag <- dt_centroids[!chain %in% c("H", "L")]
    
    # Costruisci pdb_id nel formato atteso da get_asymmetric_potential
    ch_h  <- dt_ab[chain == "H", unique(chain)]
    ch_l  <- dt_ab[chain == "L", unique(chain)]
    ch_ag <- dt_ag[, unique(chain)]
    pdb_id_str <- paste(tools::file_path_sans_ext(file_name), ch_h, ch_l, ch_ag, sep = "_")
    
    dt_ab[, .dummy := 1L]
    dt_ag[, .dummy := 1L]
    
    dt_contacts <- dt_ab[dt_ag, on = ".dummy", allow.cartesian = TRUE] |>
      _[, dist := sqrt((x - i.x)^2 + (y - i.y)^2 + (z - i.z)^2)] |>
      _[dist <= DistCutoff] |>
      _[resid %in% aa & i.resid %in% aa] |>
      _[, .(
        pdb_id   = pdb_id_str,
        resid_ab = resid,
        resno_ab = resno,
        chain_ab = chain,
        resid_ag = i.resid,
        resno_ag = i.resno,
        chain_ag = i.chain
      )]
    return(list(ok = TRUE, contacts = dt_contacts))
  }, error = function(e) {
    return(list(ok = FALSE, filename = file_name, path = pdb_path, error = e$message))
  })
}

### Esegui in parallelo
with_progress({
  results_list <- future_map(
    all_pdbs,
    gen_df_contacts,   # <- nome corretto della funzione
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
  )
})

### Filtra risultati validi
valid_results <- keep(results_list, ~ .x$ok)
failed_files  <- map_chr(discard(results_list, ~ .x$ok), "filename")
cat("File falliti:", length(failed_files), "\n")
if (length(failed_files) > 0) {
  writeLines(failed_files, file.path(results_dir, "failed_files.txt"))
}

### Combina tutti i contatti
df_contacts <- rbindlist(
  map(valid_results, "contacts"),
  use.names = TRUE
)

### Conteggi residui all'interfaccia
residue_counts_ab <- df_contacts[, .(resid_ab, resno_ab, chain_ab, pdb_id)] |>
  unique() |>
  _[, .N, by = resid_ab]

residue_counts_ag <- df_contacts[, .(resid_ag, resno_ag, chain_ag, pdb_id)] |>
  unique() |>
  _[, .N, by = resid_ag]

### Salvataggio
saveRDS(df_contacts,      file.path(results_dir, "df_contacts.rds"))
saveRDS(residue_counts_ab, file.path(results_dir, "residue_counts_ab.rds"))
saveRDS(residue_counts_ag, file.path(results_dir, "residue_counts_ag.rds"))

fwrite(df_contacts,       file.path(results_dir, "df_contacts.csv"))
fwrite(residue_counts_ab, file.path(results_dir, "residue_counts_ab.csv"))
fwrite(residue_counts_ag, file.path(results_dir, "residue_counts_ag.csv"))

### Compute frequencies and statistical potentials

# # ASYMMETRIC PART
# F(x,y) -> df_fxy
# F(x) -> df_fx
# F(y) -> df_fy
# E(x,y) = - ln(F(x,y) / (F(x) * F(y))) * 2.479
# 
# H1: 26:32
# H2: 52:56
# H3: 95:102
# 
# L1: 24:34
# L2: 50:56
# L3: 89:97
get_asymmetric_potential <- function(df_contacts, part = 'all', sigma = S) {
  aa <- aa.table$aa3[1:20]
  kBT <- 2.479
  
  all_pairs <- CJ(resid_ab = aa, resid_ag = aa) |>
    _[, pair := paste(resid_ab, resid_ag, sep = "-")] |>
    _[, .(pair)]
  
  df_contacts <- df_contacts |>
    _[(resid_ab %in% aa) & (resid_ag %in% aa)] |>
    _[, c("pdb_id", "ch_h", "ch_l", "ch_ag") := tstrsplit(pdb_id, "_")]
  
  if (part == 'all') {
    df_contacts <- df_contacts
  } else if (part == 'l1') {
    df_contacts <- df_contacts |>
      _[(resno_ab %in% c(24:34)) & (chain_ab == ch_l)]
  } else if (part == 'l2') {
    df_contacts <- df_contacts |>
      _[(resno_ab %in% c(50:56)) & (chain_ab == ch_l)]
  } else if (part == 'l3') {
    df_contacts <- df_contacts |>
      _[(resno_ab %in% c(89:97)) & (chain_ab == ch_l)]
  } else if (part == 'h1') {
    df_contacts <- df_contacts |>
      _[(resno_ab %in% c(26:32)) & (chain_ab == ch_h)]
  } else if (part == 'h2') {
    df_contacts <- df_contacts |>
      _[(resno_ab %in% c(52:56)) & (chain_ab == ch_h)]
  } else if (part == 'h3') {
    df_contacts <- df_contacts |>
      _[(resno_ab %in% c(95:102)) & (chain_ab == ch_h)]
  } else {
    stop("Invalid value for 'part'. Use one of: 'l1', 'l2', 'l3', 'h1', 'h2', 'h3', or 'all'.")
  }
  
  df_fxy <- df_contacts |>
    _[, .(count = .N), by = .(resid_ag, resid_ab)] |>
    _[, freq := count / sum(count)]
  
  df_fx <- df_contacts |>
    _[, .(count = .N), by = resid_ag] |>
    _[, freq := count / sum(count)]
  
  df_fy <- df_contacts |>
    _[, .(count = .N), by = resid_ab] |>
    _[, freq := count / sum(count)]
  
  df_potential <- df_fxy |>
    _[, pair := paste(resid_ab, resid_ag, sep = "-")] |>
    _[, .(pair, resid_ag, resid_ab, freq, count)] |>          # <- mantieni 'count'
    _[df_fx[, .(freq_ag = freq), by = resid_ag], on = "resid_ag"] |>
    _[df_fy[, .(freq_ab = freq), by = resid_ab], on = "resid_ab"] |>
    _[, potential := kBT * log(1 + count * sigma) -
        kBT * log(1 + count * sigma * (freq / (freq_ag * freq_ab)))]
  
  df_potential_complete <- merge(
    all_pairs,
    df_potential,
    by = "pair",
    all.x = TRUE
  )[, .(resid_ag, resid_ab, potential)] |>
    _[, part := part]
  
  return(df_potential_complete)
}

parts <- c("all", "h1", "h2", "h3", "l1", "l2", "l3")

potentials_list <- map(parts, ~ get_asymmetric_potential(df_contacts, part = .x))
names(potentials_list) <- parts

# Unisci tutto in un'unica tabella lunga, già con la colonna `part` per distinguerle
df_potential_combined <- rbindlist(potentials_list)

df_potential_combined

saveRDS(df_potential_combined, file.path(results_dir, "whole_int_asym_potential.rds"))

fwrite(df_potential_combined,       file.path(results_dir, "whole_int_asym_potential.csv"))


potenziale_sippl <- read.csv("/Users/lorenzosisti/Downloads/potenziali_statistici_whole_26_06_data_table_sippl/whole_int_asym_potential.csv")
potenziale_non_sippl <- read.csv("/Users/lorenzosisti/Downloads/potenziali_statistici_whole_26_06_data_table/whole_int_asym_potential.csv")

sum(is.na(potenziale_sippl$potential))
sum(is.na(potenziale_non_sippl$potential))

# Le coppie con NA sono esattamente le stesse nei due file?
na_sippl     <- which(is.na(potenziale_sippl$potential))
na_non_sippl <- which(is.na(potenziale_non_sippl$potential))
identical(na_sippl, na_non_sippl)

correlazione_potenziali <- cor(potenziale_sippl$potential, potenziale_non_sippl$potential)

# ASYMMETRIC approach: Compute amino acid frequencies for antibody (paratope) and antigen (epitope)
paratope_freq <- residue_counts_interface_ab_sum / sum(residue_counts_interface_ab_sum)
epitope_freq  <- residue_counts_interface_ligando_sum / sum(residue_counts_interface_ligando_sum)

# SYMMETRIC approach: Combine counts from both chains and compute normalized frequencies
par_plus_epi <- residue_counts_interface_ab_sum + residue_counts_interface_ligando_sum
residue_freq <- par_plus_epi / sum(par_plus_epi)

# -------------------------
# ASYMMETRIC POTENTIAL
# -------------------------
contact_freq_asym <- contact_matrix_asym_sum / sum(contact_matrix_asym_sum)

# sum(contact_matrix_asym_sum) = 169073

expected_asym <- outer(paratope_freq, epitope_freq, "*")
ratio_asym <- contact_freq_asym / expected_asym

# Computing potentials according to Sippl corrections
V_asym <- (log(1 + contact_matrix_asym_sum * S) - log(1 + contact_matrix_asym_sum * S * ratio_asym)) * 2.479

# -------------------------
# SYMMETRIC: ADJUSTING REFERENCE STATE ACCORDING TO RANDOM MIXING IN QUASI-CHEMICAL APPROXIMATION 
# -------------------------

# Keep only the lower triangle (unique unordered pairs)
mask_sym <- lower.tri(contact_matrix_sym_sum, diag = TRUE)
contact_matrix_sym_sum[!mask_sym] <- 0

# Observed frequencies on the 210 unique pairs
contact_freq_sym <- matrix(0, nrow = 20, ncol = 20,
                           dimnames = list(amino_acids, amino_acids))
contact_freq_sym[mask_sym] <- contact_matrix_sym_sum[mask_sym] / sum(contact_matrix_sym_sum[mask_sym])

# sum(contact_matrix_sym_sum[mask_sym]) = 169073

##########################################################################################################
##########################################################################################################
##########################################################################################################

# Expected frequencies for unordered pairs: MIYAZAWA-JERNIGAN 1985 IN THE QUASI-CHEMICAL APPROX I HAVE RANDOM MIXING AS THE REFERENCE
# diagonal   -> p_i^2
# off-diagonal -> 2 * p_i * p_j 
expected_sym <- outer(residue_freq, residue_freq, "*")
expected_sym[row(expected_sym) != col(expected_sym)] <- 2 * expected_sym[row(expected_sym) != col(expected_sym)] 
expected_sym[!mask_sym] <- 0

# Observed / expected ratio
ratio_sym <- matrix(0, nrow = 20, ncol = 20,
                    dimnames = list(amino_acids, amino_acids))
ratio_sym[mask_sym] <- contact_freq_sym[mask_sym] / expected_sym[mask_sym]

# Sparse-data Sippl potential
V_sym <- matrix(0, nrow = 20, ncol = 20,
                dimnames = list(amino_acids, amino_acids))
V_sym[mask_sym] <- (log(1 + contact_matrix_sym_sum[mask_sym] * S) - log(1 + contact_matrix_sym_sum[mask_sym] * S * ratio_sym[mask_sym])) * 2.479

### Prepare and plot heatmaps
# Compute global min/max across both matrices for consistent color scaling
pot_min <- abs(min(min(V_sym), min(V_asym)))
pot_max <- abs(max(max(V_sym), max(V_asym)))
# Assign row/column names for the asymmetric heatmap
rownames(V_asym) <- paste(amino_acids, "Ab", sep = "_")
colnames(V_asym) <- paste(amino_acids, "Ag", sep = "_")

p_asym <- pheatmap(V_asym,
                   color = colorRampPalette(c("gold1", "white", "dodgerblue2"))(50),
                   breaks = seq(- max(pot_min, pot_max), max(pot_min, pot_max), length.out = 51),
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   main = paste("Asymmetric Potential"),
                   display_numbers = FALSE,
                   fontsize = 10)
#silent = TRUE) 
heatmap_asym_plots <- as.ggplot(p_asym$gtable)

p_sym <- pheatmap(V_sym,
                  color = colorRampPalette(c("gold1", "white", "dodgerblue2"))(50),
                  breaks = seq(- max(pot_min, pot_max), max(pot_min, pot_max), length.out = 51),
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  main = paste("Symmetric Potential"),
                  display_numbers = FALSE,
                  fontsize = 10)
#silent = TRUE) 
heatmap_sym_plots <- as.ggplot(p_sym$gtable) 

# Save as .csv
write.csv(V_asym, file = file.path(results_dir, paste0("V_asym.csv")))
write.csv(V_sym,  file = file.path(results_dir, paste0("V_sym.csv")))

# Prepare and export statistical potentials data frames for further analysis

sym_df <- melt(V_sym)
idx <- which(lower.tri(V_sym, diag = TRUE), arr.ind = TRUE)
aa1 <- rownames(V_sym)[idx[,1]]
aa2 <- colnames(V_sym)[idx[,2]]
pair <- paste(pmin(aa1, aa2), pmax(aa1, aa2), sep = "-")
val  <- V_sym[idx]
sym_pairs <- data.frame(pair = pair, value = as.numeric(val), stringsAsFactors = FALSE)
sym_pairs_sep <- sym_pairs %>%
  separate(pair, into = c("aa1", "aa2"), sep = "-")

write.csv(sym_pairs_sep, file = file.path(results_dir, paste0("V_sym_df.csv")))


asym_df <- melt(V_asym) %>%
  rename(
    aa1 = Var1,
    aa2 = Var2
  )

write.csv(asym_df, file = file.path(results_dir, paste0("V_asym_df.csv")))

# Combine both heatmaps into a single plot
(heatmap_asym_plots | heatmap_sym_plots)

### Check cose di Leonardo nuovo dataset EM ###

sym_df <- melt(V_sym)
leonardo_V_sym <- read.csv("/Users/lorenzosisti/Downloads/df_symmetric_statistical_potentials_centroids_minimized.csv")

# Filtra e tieni solo le colonne di interesse
df_all <- leonardo_V_sym %>%
  filter(part == "all") %>%
  select(pair, potential)

head(df_all)

## 1) df_all: hai già fatto filter(part=="all") e select(pair, potential)
##    Assicurati che 'pair' sia tipo "ALA-ARG" (ordine alfabetico con '-')

## 2) Estrai le 210 coppie da V_sym
stopifnot(is.matrix(V_sym), !is.null(rownames(V_sym)), !is.null(colnames(V_sym)))

idx <- which(lower.tri(V_sym, diag = TRUE), arr.ind = TRUE)
aa1 <- rownames(V_sym)[idx[,1]]
aa2 <- colnames(V_sym)[idx[,2]]

## Ordina alfabeticamente ogni coppia per avere una chiave canonica "AAA-BBB"
pair <- paste(pmin(aa1, aa2), pmax(aa1, aa2), sep = "-")
val  <- V_sym[idx]

sym_pairs <- data.frame(pair = pair, value = as.numeric(val), stringsAsFactors = FALSE)

sym_pairs_sep <- sym_pairs %>%
  separate(pair, into = c("aa1", "aa2"), sep = "-")

head(sym_pairs_sep)
#write.csv(sym_pairs_sep, file = file.path(results_dir, paste0("V_sym_df.csv")))


## 3) Unisci con il CSV filtrato
joined <- merge(df_all, sym_pairs, by = "pair", all = FALSE)

## 4) Correlazioni
cor_pearson  <- cor(joined$potential, joined$value, method = "pearson",  use = "complete.obs")

###

asym_df <- melt(V_asym)

#write.csv(asym_df, file = file.path(results_dir, paste0("V_asym_df.csv")))
leonardo_V_asym <- read.csv("/Users/lorenzosisti/Downloads/df_asymmetric_statistical_potentials_centroids_minimized.csv")
leonardo_V_asym$Var1 <- paste(leonardo_V_asym$resid_ab, "_Ab", sep = "")
leonardo_V_asym$Var2 <- paste(leonardo_V_asym$resid_ag, "_Ag", sep = "")
leonardo_V_asym <- leonardo_V_asym[leonardo_V_asym$part == "all",]


# Assumendo che asym_df sia già creato con melt(V_asym)
# E che leonardo_V_asym abbia già le colonne Var1 e Var2

# Merge sui nomi delle coppie ordinate
merged_df <- merge(asym_df, 
                   leonardo_V_asym[, c("Var1", "Var2", "potential")],
                   by = c("Var1", "Var2"))

# Calcolo della correlazione di Pearson
correlation <- cor(merged_df$value, merged_df$potential, method = "pearson")

(heatmap_asym_plots | heatmap_sym_plots)



