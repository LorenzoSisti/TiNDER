#########################################################################################
# THIS SCRIPT COMPUTES AMINO ACID CONTACT MATRICES AND WHOLE-INTERFACE STATISTICAL POTENTIALS FROM A DIRECTORY OF PDB STRUCTURES. 
# IT ANALYZES ANTIBODY–ANTIGEN INTERACTIONS BY:
# 1) CALCULATING SIDE-CHAIN CENTROIDS AND CONTACT MATRICES
# 2) IDENTIFYING INTERFACE RESIDUES BETWEEN ANTIBODY AND ANTIGEN
# 3) SUMMARIZING CONTACT FREQUENCIES ACROSS ALL STRUCTURES
# 4) COMPUTING ASYMMETRIC AND SYMMETRIC WHOLE-INTERFACE STATISTICAL POTENTIALS
# 5) VISUALIZING THE RESULTS AS HEATMAPS.
#########################################################################################

### Required libraries
pacman::p_load(bio3d, dplyr, future, furrr, purrr, progressr,
               pheatmap, patchwork, ggplotify, reshape2, tidyr, data.table, ggplot2)

# Define the path to a custom function files
#source("/path/to/your/custom/functions/files/functions.R")
source("/Users/lorenzosisti/Documents/Script_ottimizzati_funzioni/functions.R")

### Set up parallelization to speed up computation
plan(multisession, workers = parallel::detectCores() - 1)
handlers(global = TRUE)
handlers("rstudio")

### Define directories and global parameters
pdb_dir <- "/Users/lorenzosisti/Downloads/database_settembre_renamed/"
results_dir <- "/Users/lorenzosisti/Downloads/potenziali_statistici_whole_29_06_data_table_sippl/"
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
  df_contacts <- copy(df_contacts)
  aa <- aa.table$aa3[1:20]
  kBT <- 2.479
  
  all_pairs <- CJ(resid_ab = aa, resid_ag = aa) |>
    _[, pair := paste(resid_ab, resid_ag, sep = "-")]
  
  df_contacts <- df_contacts |>
    _[(resid_ab %in% aa) & (resid_ag %in% aa)] |>
    _[, c("pdb_id", "ch_h", "ch_l", "ch_ag") := tstrsplit(pdb_id, "_")]
  
  if (part == 'all') {
    df_contacts <- df_contacts
  } else if (part == 'l1') {
    df_contacts <- df_contacts[(resno_ab %in% c(24:34)) & (chain_ab == ch_l)]
  } else if (part == 'l2') {
    df_contacts <- df_contacts[(resno_ab %in% c(50:56)) & (chain_ab == ch_l)]
  } else if (part == 'l3') {
    df_contacts <- df_contacts[(resno_ab %in% c(89:97)) & (chain_ab == ch_l)]
  } else if (part == 'h1') {
    df_contacts <- df_contacts[(resno_ab %in% c(26:32)) & (chain_ab == ch_h)]
  } else if (part == 'h2') {
    df_contacts <- df_contacts[(resno_ab %in% c(52:56)) & (chain_ab == ch_h)]
  } else if (part == 'h3') {
    df_contacts <- df_contacts[(resno_ab %in% c(95:102)) & (chain_ab == ch_h)]
  } else {
    stop("Invalid value for 'part'. Use one of: 'l1', 'l2', 'l3', 'h1', 'h2', 'h3', or 'all'.")
  }
  
  df_fxy <- df_contacts[, .(count = .N), by = .(resid_ag, resid_ab)][, freq := count / sum(count)]
  df_fx  <- df_contacts[, .(count = .N), by = resid_ag][, freq := count / sum(count)]
  df_fy  <- df_contacts[, .(count = .N), by = resid_ab][, freq := count / sum(count)]
  
  df_potential <- df_fxy |>
    _[, pair := paste(resid_ab, resid_ag, sep = "-")] |>
    _[, .(pair, resid_ag, resid_ab, freq, count)] |>
    _[df_fx[, .(freq_ag = freq), by = resid_ag], on = "resid_ag"] |>
    _[df_fy[, .(freq_ab = freq), by = resid_ab], on = "resid_ab"] |>
    _[, potential := kBT * log(1 + count * sigma) -
        kBT * log(1 + count * sigma * (freq / (freq_ag * freq_ab)))]
  
  df_potential_complete <- merge(
    all_pairs,
    df_potential[, .(pair, potential)],
    by = "pair",
    all.x = TRUE
  )[, .(resid_ag, resid_ab, potential)] |>
    _[, potential := fifelse(is.na(potential), 0, potential)] |>
    _[, part := part]
  
  return(df_potential_complete)
}

get_symmetric_potential <- function(df_contacts, part = 'all', sigma = S) {
  df_contacts <- copy(df_contacts)
  aa <- aa.table$aa3[1:20]
  kBT <- 2.479
  
  # Tutte le 210 coppie non ordinate (i <= j alfabeticamente, evita duplicati ARG-LYS / LYS-ARG)
  all_pairs <- CJ(resid_i = aa, resid_j = aa) |>
    _[resid_i <= resid_j] |>
    _[, pair := paste(resid_i, resid_j, sep = "-")]
  
  df_contacts <- df_contacts |>
    _[(resid_ab %in% aa) & (resid_ag %in% aa)] |>
    _[, c("pdb_id", "ch_h", "ch_l", "ch_ag") := tstrsplit(pdb_id, "_")]
  
  if (part == 'all') {
    df_contacts <- df_contacts
  } else if (part == 'l1') {
    df_contacts <- df_contacts[(resno_ab %in% c(24:34)) & (chain_ab == ch_l)]
  } else if (part == 'l2') {
    df_contacts <- df_contacts[(resno_ab %in% c(50:56)) & (chain_ab == ch_l)]
  } else if (part == 'l3') {
    df_contacts <- df_contacts[(resno_ab %in% c(89:97)) & (chain_ab == ch_l)]
  } else if (part == 'h1') {
    df_contacts <- df_contacts[(resno_ab %in% c(26:32)) & (chain_ab == ch_h)]
  } else if (part == 'h2') {
    df_contacts <- df_contacts[(resno_ab %in% c(52:56)) & (chain_ab == ch_h)]
  } else if (part == 'h3') {
    df_contacts <- df_contacts[(resno_ab %in% c(95:102)) & (chain_ab == ch_h)]
  } else {
    stop("Invalid value for 'part'. Use one of: 'l1', 'l2', 'l3', 'h1', 'h2', 'h3', or 'all'.")
  }
  
  # Etichetta non ordinata per ogni contatto: resid_i <= resid_j alfabeticamente
  df_contacts <- df_contacts |>
    _[, `:=`(
      resid_i = fifelse(resid_ab <= resid_ag, resid_ab, resid_ag),
      resid_j = fifelse(resid_ab <= resid_ag, resid_ag, resid_ab)
    )]
  
  # F(x,y): frequenza osservata della coppia non ordinata
  df_fxy <- df_contacts[, .(count = .N), by = .(resid_i, resid_j)][, freq := count / sum(count)]
  
  # F(x): UNICA distribuzione marginale, ottenuta mettendo insieme ab e ag
  df_pool <- rbindlist(list(
    df_contacts[, .(resid = resid_ab)],
    df_contacts[, .(resid = resid_ag)]
  ))
  df_f <- df_pool[, .(count = .N), by = resid][, freq := count / sum(count)]
  
  freq_i_dt <- df_f[, .(resid_i = resid, freq_i = freq)]
  freq_j_dt <- df_f[, .(resid_j = resid, freq_j = freq)]
  
  df_potential <- df_fxy |>
    _[, pair := paste(resid_i, resid_j, sep = "-")] |>
    _[, .(pair, resid_i, resid_j, freq, count)] |>
    _[freq_i_dt, on = "resid_i"] |>
    _[freq_j_dt, on = "resid_j"] |>
    _[, ref_freq := fifelse(resid_i == resid_j, freq_i * freq_j, 2 * freq_i * freq_j)] |>
    _[, potential := kBT * log(1 + count * sigma) -
        kBT * log(1 + count * sigma * (freq / ref_freq))]
  
  df_potential_complete <- merge(
    all_pairs,
    df_potential[, .(pair, potential)],
    by = "pair",
    all.x = TRUE
  )[, .(resid_i, resid_j, potential)] |>
    _[, potential := fifelse(is.na(potential), 0, potential)] |>
    _[, part := part]
  
  return(df_potential_complete)
}

parts <- c("all", "h1", "h2", "h3", "l1", "l2", "l3")

### ASYM ###

potentials_list <- map(parts, ~ get_asymmetric_potential(df_contacts, part = .x))
names(potentials_list) <- parts

# Unisci tutto in un'unica tabella lunga, già con la colonna `part` per distinguerle
df_potential_combined <- rbindlist(potentials_list)

df_potential_combined

saveRDS(df_potential_combined, file.path(results_dir, "whole_int_asym_potential.rds"))

fwrite(df_potential_combined,       file.path(results_dir, "whole_int_asym_potential.csv"))

### SYM ###

sym_potentials_list <- map(parts, ~ get_symmetric_potential(df_contacts, part = .x))
names(sym_potentials_list) <- parts

# Unisci tutto in un'unica tabella lunga, già con la colonna `part` per distinguerle
df_sym_potential_combined <- rbindlist(sym_potentials_list)

df_sym_potential_combined

saveRDS(df_sym_potential_combined, file.path(results_dir, "whole_int_sym_potential.rds"))

fwrite(df_sym_potential_combined,       file.path(results_dir, "whole_int_sym_potential.csv"))


### Costruisce la matrice 20x20 a partire dal data.table long
build_potential_matrix <- function(df_potential, part_name, symmetric = FALSE, aa_order = amino_acids) {
  
  dt <- df_potential[part == part_name]
  
  if (!symmetric) {
    # Caso asimmetrico: righe = resid_ab, colonne = resid_ag
    mat_dt <- dcast(dt, resid_ab ~ resid_ag, value.var = "potential")
    rn  <- mat_dt$resid_ab
    mat <- as.matrix(mat_dt[, -1, with = FALSE])
    rownames(mat) <- rn
  } else {
    # Caso simmetrico: la tabella ha solo resid_i <= resid_j -> va "specchiata"
    dt_full <- rbindlist(list(
      dt[, .(resid_i, resid_j, potential)],
      dt[resid_i != resid_j, .(resid_i = resid_j, resid_j = resid_i, potential)]
    ))
    mat_dt <- dcast(dt_full, resid_i ~ resid_j, value.var = "potential")
    rn  <- mat_dt$resid_i
    mat <- as.matrix(mat_dt[, -1, with = FALSE])
    rownames(mat) <- rn
  }
  
  # Riordina righe/colonne secondo la scala Kyte-Doolittle
  mat <- mat[aa_order, aa_order]
  
  return(mat)
}

### Genera (e opzionalmente salva) l'heatmap per un singolo 'part'
plot_potential_heatmap <- function(df_potential,
                                   part_name,
                                   symmetric  = FALSE,
                                   aa_order   = amino_acids,
                                   title_prefix = "Whole-Interface Potential",
                                   save_dir   = results_dir,
                                   save       = TRUE,
                                   width = 7, height = 6) {
  
  mat <- build_potential_matrix(df_potential, part_name, symmetric = symmetric, aa_order = aa_order)
  
  pot_min <- min(mat, na.rm = TRUE)
  pot_max <- max(mat, na.rm = TRUE)
  lim     <- max(abs(pot_min), abs(pot_max))
  
  type_label <- if (symmetric) "Symmetric" else "Asymmetric"
  main_title <- paste(type_label, title_prefix, "-", toupper(part_name))
  
  p <- pheatmap(mat,
                color = colorRampPalette(c("gold1", "white", "dodgerblue2"))(50),
                breaks = seq(-lim, lim, length.out = 51),
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                main = main_title,
                display_numbers = FALSE,
                fontsize = 10,
                silent = TRUE)
  
  gg_heatmap <- as.ggplot(p$gtable)
  
  if (save) {
    file_suffix <- if (symmetric) "sym" else "asym"
    file_name   <- paste0("heatmap_", file_suffix, "_", part_name, ".png")
    ggsave(filename = file.path(save_dir, file_name),
           plot = gg_heatmap, width = width, height = height, dpi = 300, bg = "white")
  }
  
  return(gg_heatmap)
}


### Heatmap per il potenziale asimmetrico
heatmaps_asym <- map(parts, ~ plot_potential_heatmap(df_potential_combined,
                                                     part_name = .x,
                                                     symmetric = FALSE))
names(heatmaps_asym) <- parts

### Heatmap per il potenziale simmetrico
heatmaps_sym <- map(parts, ~ plot_potential_heatmap(df_sym_potential_combined,
                                                    part_name = .x,
                                                    symmetric = TRUE))
names(heatmaps_sym) <- parts
