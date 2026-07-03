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
results_dir <- "/Users/lorenzosisti/Downloads/potenziali_statistici_gr_30_06_data_table_sippl/"
dir.create(results_dir, showWarnings = FALSE)

# Distance cutoff (Å) to define contact between side-chains centroids
DistCutoff <- 8.5  
S <- 0.02
Nsteps <- 3

# Amino acids ordered according to the Kyte-Doolittle hydrophobicity scale
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE") 

set.seed(1234)

all_pdbs <- list.files(pdb_dir, pattern = "*.pdb", recursive = TRUE, full.names = TRUE)

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
    
    ch_h  <- dt_ab[chain == "H", unique(chain)]
    ch_l  <- dt_ab[chain == "L", unique(chain)]
    ch_ag <- dt_ag[, unique(chain)]
    pdb_id_str <- paste(tools::file_path_sans_ext(file_name), ch_h, ch_l, ch_ag, sep = "_")
    
    dt_ab[, .dummy := 1L]
    dt_ag[, .dummy := 1L]
    
    # Contatti grezzi, con le coordinate dei due residui (servono per il centroide del BS)
    dt_contacts <- dt_ab[dt_ag, on = ".dummy", allow.cartesian = TRUE] |>
      _[, dist := sqrt((x - i.x)^2 + (y - i.y)^2 + (z - i.z)^2)] |>
      _[dist <= DistCutoff] |>
      _[resid %in% aa & i.resid %in% aa] |>
      _[, .(
        pdb_id   = pdb_id_str,
        resid_ab = resid,  resno_ab = resno,  chain_ab = chain,
        x_ab = x, y_ab = y, z_ab = z,
        resid_ag = i.resid, resno_ag = i.resno, chain_ag = i.chain,
        x_ag = i.x, y_ag = i.y, z_ag = i.z
      )]
    
    if (nrow(dt_contacts) == 0) {
      return(list(ok = FALSE, filename = file_name, path = pdb_path, error = "no contacts found"))
    }
    
    ### --- ASSEGNAZIONE AI RING (centroide del binding site) ---
    
    # Residui di interfaccia unici, lato anticorpo e lato antigene
    dt_if_ab <- unique(dt_contacts[, .(chain = chain_ab, resno = resno_ab, resid = resid_ab, x = x_ab, y = y_ab, z = z_ab)])
    dt_if_ag <- unique(dt_contacts[, .(chain = chain_ag, resno = resno_ag, resid = resid_ag, x = x_ag, y = y_ag, z = z_ag)])
    dt_if_all <- rbindlist(list(dt_if_ab, dt_if_ag))
    
    center_BS <- dt_if_all[, .(x = mean(x), y = mean(y), z = mean(z))]
    
    dt_if_all[, dist_centroid := sqrt((x - center_BS$x)^2 + (y - center_BS$y)^2 + (z - center_BS$z)^2)]
    
    r_max   <- max(dt_if_all$dist_centroid)
    borders <- c(0, 0.37, 0.64, 1) * r_max
    
    dt_if_all[, ring := findInterval(dist_centroid, borders, left.open = TRUE)]
    
    ring_lookup <- dt_if_all[, .(chain, resno, ring)]
    
    # Assegna il ring a ciascun lato del contatto
    dt_contacts <- merge(dt_contacts, ring_lookup,
                         by.x = c("chain_ab", "resno_ab"), by.y = c("chain", "resno"), all.x = TRUE)
    setnames(dt_contacts, "ring", "ring_ab")
    
    dt_contacts <- merge(dt_contacts, ring_lookup,
                         by.x = c("chain_ag", "resno_ag"), by.y = c("chain", "resno"), all.x = TRUE)
    setnames(dt_contacts, "ring", "ring_ag")
    
    # Scarta eventuali contatti con ring fuori range (non dovrebbe succedere, ma per sicurezza)
    # dt_contacts <- dt_contacts[ring_ab %in% 1:Nsteps & ring_ag %in% 1:Nsteps]
    
    # Etichetta non ordinata della coppia di ring: (1,2) e (2,1) -> "1-2"
    dt_contacts[, `:=`(
      ring_i = pmin(ring_ab, ring_ag),
      ring_j = pmax(ring_ab, ring_ag)
    )]
    dt_contacts[, ring_pair := paste(ring_i, ring_j, sep = "-")]
    
    dt_contacts <- dt_contacts[, .(pdb_id, resid_ab, resno_ab, chain_ab,
                                   resid_ag, resno_ag, chain_ag,
                                   ring_ab, ring_ag, ring_pair)]
    
    return(list(ok = TRUE, contacts = dt_contacts))
  }, error = function(e) {
    return(list(ok = FALSE, filename = file_name, path = pdb_path, error = e$message))
  })
}

### Esegui in parallelo
with_progress({
  results_list <- future_map(
    all_pdbs,
    gen_df_contacts,
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

### Combina tutti i contatti (tabella unica, con colonna ring_pair)
df_contacts <- rbindlist(map(valid_results, "contacts"), use.names = TRUE)

saveRDS(df_contacts, file.path(results_dir, "df_contacts_all_rings.rds"))
fwrite(df_contacts,  file.path(results_dir, "df_contacts_all_rings.csv"))

### --- SPLIT IN 6 DATAFRAME, UNO PER RING PAIR ---
ring_pairs <- c("1-1", "2-2", "3-3", "1-2", "2-3", "1-3")

df_contacts_by_ring <- map(ring_pairs, ~ df_contacts[ring_pair == .x])
names(df_contacts_by_ring) <- ring_pairs

# Salvataggio: sia come lista RDS unica, sia come singoli RDS/CSV
saveRDS(df_contacts_by_ring, file.path(results_dir, "df_contacts_by_ring_list.rds"))

for (pair in ring_pairs) {
  saveRDS(df_contacts_by_ring[[pair]],
          file.path(results_dir, paste0("df_contacts_ring_", pair, ".rds")))
  fwrite(df_contacts_by_ring[[pair]],
         file.path(results_dir, paste0("df_contacts_ring_", pair, ".csv")))
}

### Conteggi residui all'interfaccia (sull'intero df_contacts, come prima)
residue_counts_ab <- df_contacts[, .(resid_ab, resno_ab, chain_ab, pdb_id)] |>
  unique() |> _[, .N, by = resid_ab]

residue_counts_ag <- df_contacts[, .(resid_ag, resno_ag, chain_ag, pdb_id)] |>
  unique() |> _[, .N, by = resid_ag]

saveRDS(residue_counts_ab, file.path(results_dir, "residue_counts_ab.rds"))
saveRDS(residue_counts_ag, file.path(results_dir, "residue_counts_ag.rds"))
fwrite(residue_counts_ab, file.path(results_dir, "residue_counts_ab.csv"))
fwrite(residue_counts_ag, file.path(results_dir, "residue_counts_ag.csv"))

get_asymmetric_potential <- function(df_contacts, ring = 'all', sigma = S) {
  df_contacts <- copy(df_contacts)
  aa <- aa.table$aa3[1:20]
  kBT <- 2.479
  
  all_pairs <- CJ(resid_ab = aa, resid_ag = aa) |>
    _[, pair := paste(resid_ab, resid_ag, sep = "-")]
  
  df_contacts <- df_contacts[(resid_ab %in% aa) & (resid_ag %in% aa)]
  
  if (ring != 'all') {
    df_contacts <- df_contacts[ring_pair == ring]
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
    all_pairs, df_potential[, .(pair, potential)], by = "pair", all.x = TRUE
  )[, .(resid_ag, resid_ab, potential)] |>
    _[, potential := fifelse(is.na(potential), 0, potential)] |>
    _[, ring_pair := ring]
  
  return(df_potential_complete)
}

get_symmetric_potential <- function(df_contacts, ring = 'all', sigma = S) {
  df_contacts <- copy(df_contacts)
  aa <- aa.table$aa3[1:20]
  kBT <- 2.479
  
  all_pairs <- CJ(resid_i = aa, resid_j = aa) |>
    _[resid_i <= resid_j] |>
    _[, pair := paste(resid_i, resid_j, sep = "-")]
  
  df_contacts <- df_contacts[(resid_ab %in% aa) & (resid_ag %in% aa)]
  
  if (ring != 'all') {
    df_contacts <- df_contacts[ring_pair == ring]
  }
  
  df_contacts <- df_contacts |>
    _[, `:=`(
      resid_i = fifelse(resid_ab <= resid_ag, resid_ab, resid_ag),
      resid_j = fifelse(resid_ab <= resid_ag, resid_ag, resid_ab)
    )]
  
  df_fxy <- df_contacts[, .(count = .N), by = .(resid_i, resid_j)][, freq := count / sum(count)]
  
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
    all_pairs, df_potential[, .(pair, potential)], by = "pair", all.x = TRUE
  )[, .(resid_i, resid_j, potential)] |>
    _[, potential := fifelse(is.na(potential), 0, potential)] |>
    _[, ring_pair := ring]
  
  return(df_potential_complete)
}

ring_pairs <- c("1-1", "2-2", "3-3", "1-2", "2-3", "1-3")

### ASYM ###
potentials_list <- map(ring_pairs, ~ get_asymmetric_potential(df_contacts, ring = .x))
names(potentials_list) <- ring_pairs
df_potential_combined <- rbindlist(potentials_list)

saveRDS(df_potential_combined, file.path(results_dir, "ring_asym_potential.rds"))
fwrite(df_potential_combined, file.path(results_dir, "ring_asym_potential.csv"))

### SYM ###
sym_potentials_list <- map(ring_pairs, ~ get_symmetric_potential(df_contacts, ring = .x))
names(sym_potentials_list) <- ring_pairs
df_sym_potential_combined <- rbindlist(sym_potentials_list)

saveRDS(df_sym_potential_combined, file.path(results_dir, "ring_sym_potential.rds"))
fwrite(df_sym_potential_combined, file.path(results_dir, "ring_sym_potential.csv"))

### Costruisce la matrice 20x20 a partire dal data.table long
build_potential_matrix <- function(df_potential, ring_name, symmetric = FALSE, aa_order = amino_acids) {
  
  dt <- df_potential[ring_pair == ring_name]
  
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

plot_potential_heatmap <- function(df_potential,
                                   ring_name,
                                   symmetric  = FALSE,
                                   aa_order   = amino_acids,
                                   title_prefix = "Radial Potential",
                                   save_dir   = results_dir,
                                   save       = TRUE,
                                   width = 7, height = 6,
                                   global_lim = NULL) {
  
  mat <- build_potential_matrix(df_potential, ring_name, symmetric = symmetric, aa_order = aa_order)
  
  if (is.null(global_lim)) {
    pot_min <- min(mat, na.rm = TRUE)
    pot_max <- max(mat, na.rm = TRUE)
    lim     <- max(abs(pot_min), abs(pot_max))
  } else {
    lim <- global_lim
  }
  
  type_label <- if (symmetric) "Symmetric" else "Asymmetric"
  main_title <- paste(type_label, title_prefix, "-", toupper(ring_name))
  
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
    file_name   <- paste0("heatmap_", file_suffix, "_", ring_name, ".png")
    ggsave(filename = file.path(save_dir, file_name),
           plot = gg_heatmap, width = width, height = height, dpi = 300, bg = "white")
  }
  
  return(gg_heatmap)
}


# Shared limit across all 6 asymmetric ring pairs
lim_asym <- max(abs(df_potential_combined$potential), na.rm = TRUE)

# Shared limit across all 6 symmetric ring pairs
lim_sym  <- max(abs(df_sym_potential_combined$potential), na.rm = TRUE)

heatmaps_asym <- map(ring_pairs, ~ plot_potential_heatmap(df_potential_combined,
                                                          ring_name = .x,
                                                          symmetric = FALSE,
                                                          global_lim = lim_asym))
names(heatmaps_asym) <- ring_pairs

heatmaps_sym <- map(ring_pairs, ~ plot_potential_heatmap(df_sym_potential_combined,
                                                         ring_name = .x,
                                                         symmetric = TRUE,
                                                         global_lim = lim_sym))
names(heatmaps_sym) <- ring_pairs
