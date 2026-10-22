### TODO: Adding a script that,based on the symmetric whole-interface calculation scripts, computes a saturation script
### I need to start from scratch as the brutta copia was definitely illeggibile ( mi sono pure perso dei pezzi vb)

#########################################################################################
# CONVERGENZA DEL POTENZIALE STATISTICO SIMMETRICO WHOLE-INTERFACE
# 1) Calcola i contatti (SC-centroid) una sola volta per tutti i PDB
# 2) Per ogni seed e per ogni numero di partizioni n, divide il dataset in n gruppi
# 3) Per ogni gruppo calcola il potenziale simmetrico whole-interface (part = "all")
#########################################################################################

### Required libraries
pacman::p_load(bio3d, dplyr, future, furrr, purrr, progressr,
               tidyr, data.table, ggplot2, latex2exp)

# Define the path to a custom function files
source("/Users/lorenzosisti/TiNDER/functions.R")

### Set up parallelization to speed up computation
plan(multisession, workers = parallel::detectCores() - 1)
handlers(global = TRUE)
handlers("rstudio")

### Define directories and global parameters
pdb_dir <- "/Users/lorenzosisti/Downloads/database_settembre_renamed/"
results_dir <- "/Users/lorenzosisti/TiNDER/data/whole_int_sym_pot_convergence"
dir.create(results_dir, showWarnings = FALSE)

S <- 0.02
repr <- "sc_centroid"
REPR_PARAMS <- list(
  sc_centroid = list(cutoff = 8.5)
)

amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE")

seed_list <- c(101, 123, 202, 303, 404, 505, 606, 707, 808, 909)
partitioning_numbers <- c(2, 3, 5, 8, 13, 21, 34, 55, 89, 144)
pdbL <- list.files(pdb_dir, pattern = "\\.pdb$", full.names = TRUE)


### Generate a contacts dataframe (identica allo script whole-interface)
### N.B. This function relies on the renumber_ab_chains() function (see functions.R)
gen_df_contacts <- function(pdb_path, representation = "sc_centroid") {
  aa <- aa.table$aa3[1:20]
  cutoff <- REPR_PARAMS[[representation]]$cutoff
  file_name <- basename(pdb_path)
  tryCatch({
    pdb_aus <- read.pdb(pdb_path)
    renumbered_df <- renumber_ab_chains(pdb_aus, pdb_path, log_file = "errors.log",
                                        representation = representation)
    if (!renumbered_df$ok) {
      return(list(ok = FALSE, filename = file_name, path = pdb_path, error = renumbered_df$error))
    }
    dt_coord <- as.data.table(renumbered_df$df_coord_renumbered)
    
    if (representation == "all_atom") {
      dt_points <- dt_coord[, .(chain, resno, resid, region, x, y, z)]
    } else {
      dt_points <- dt_coord[, .(x = mean(x), y = mean(y), z = mean(z)),
                            by = .(chain, resno, resid, region)]
    }
    dt_ab <- dt_points[chain %in% c("H", "L")]
    dt_ag <- dt_points[!chain %in% c("H", "L")]
    
    # Build pdb_id
    ch_h  <- dt_ab[chain == "H", unique(chain)]
    ch_l  <- dt_ab[chain == "L", unique(chain)]
    ch_ag <- dt_ag[, unique(chain)]
    pdb_id_str <- paste(tools::file_path_sans_ext(file_name), ch_h, ch_l, ch_ag, sep = "_")
    
    dt_contacts <- gen_unique_residue_contacts(dt_ab, dt_ag, cutoff)
    if (is.null(dt_contacts)) {
      return(list(ok = FALSE, filename = file_name, path = pdb_path, error = "no contacts"))
    }
    dt_contacts <- dt_contacts[resid_ab %in% aa & resid_ag %in% aa] |>
      _[, pdb_id := pdb_id_str] |>
      setcolorder(c("pdb_id", "resid_ab", "resno_ab", "chain_ab", "region_ab",
                    "resid_ag", "resno_ag", "chain_ag"))
    
    return(list(ok = TRUE, contacts = dt_contacts))
  }, error = function(e) {
    return(list(ok = FALSE, filename = file_name, path = pdb_path, error = e$message))
  })
}


### Symmetric potential (identica allo script whole-interface)
get_symmetric_potential <- function(df_contacts, part = 'all', sigma = S) {
  df_contacts <- copy(df_contacts)
  aa <- aa.table$aa3[1:20]
  kBT <- 2.479
  
  # All the 210 non ordered amino acids pairs (avoid ARG-LYS / LYS-ARG duplicates)
  all_pairs <- CJ(resid_i = aa, resid_j = aa) |>
    _[resid_i <= resid_j] |>
    _[, pair := paste(resid_i, resid_j, sep = "-")]
  
  df_contacts <- df_contacts[(resid_ab %in% aa) & (resid_ag %in% aa)]
  
  if (part == 'all') {
    df_contacts <- df_contacts
  } else {
    df_contacts <- df_contacts[region_ab == part]
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


################################################################################
### 1) CONTATTI PER TUTTI I PDB (calcolati una sola volta)
################################################################################

cat("\n>>> Representation:", repr, "| cutoff:", REPR_PARAMS[[repr]]$cutoff, "A\n")

with_progress({
  results_list <- future_map(
    pdbL,
    gen_df_contacts,
    representation = repr,
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
  )
})
names(results_list) <- pdbL

valid_results <- keep(results_list, ~ .x$ok)
failed_files  <- map_chr(discard(results_list, ~ .x$ok), "filename")
cat("File falliti:", length(failed_files), "\n")
if (length(failed_files) > 0) {
  writeLines(failed_files, file.path(results_dir, "failed_files.txt"))
}

# idcol: tiene traccia del percorso del PDB di origine, serve per filtrare per gruppo
df_contacts <- rbindlist(
  map(valid_results, "contacts"),
  use.names = TRUE,
  idcol = "pdb_path"
)
fwrite(df_contacts, file.path(results_dir, "df_contacts.csv"))
setkey(df_contacts, pdb_path)

# Usiamo solo i PDB validi per costruire i gruppi
pdbL <- names(valid_results)


################################################################################
### 2) SPLIT DEL DATASET + POTENZIALE SIMMETRICO PER OGNI GRUPPO
################################################################################

for (i in seq_along(seed_list)) {
  
  current_seed <- seed_list[i]
  
  for (n in partitioning_numbers) {
    
    set.seed(current_seed)
    shuffled_files <- sample(pdbL)
    
    # Calcola quanti file per gruppo (floor)
    files_per_group <- floor(length(shuffled_files) / n)
    
    # Prendi solo i file che useremo (gli altri verranno ignorati)
    usable_files <- shuffled_files[1:(files_per_group * n)]
    
    # Crea i gruppi (lista di vettori di percorsi)
    split_groups <- split(usable_files, rep(1:n, each = files_per_group))
    
    # Ogni run con seed diversa avrà una directory dedicata
    split_dir <- file.path(results_dir, paste0("seed_", current_seed), paste0("split_", n))
    dir.create(split_dir, recursive = TRUE, showWarnings = FALSE)
    
    sym_potentials_split <- vector("list", length(split_groups))
    
    for (j in seq_along(split_groups)) {
      
      # 1. Scrivi l'elenco dei PDB del gruppo (un percorso per riga)
      files_for_this_group <- split_groups[[j]]
      writeLines(files_for_this_group,
                 con = file.path(split_dir, paste0("group_", j, ".txt")))
      
      # 2. Contatti dei soli PDB di questo gruppo
      df_contacts_group <- df_contacts[pdb_path %in% files_for_this_group]
      
      # 3. Potenziale simmetrico whole-interface del gruppo
      sym_potentials_split[[j]] <- get_symmetric_potential(df_contacts_group, part = "all") |>
        _[, `:=`(seed = current_seed, split = n, group = j)]
    }
    
    # Potenziali di tutti i gruppi di questo split in un unico dataframe long
    df_sym_potential_split <- rbindlist(sym_potentials_split)
    saveRDS(df_sym_potential_split, file.path(split_dir, "whole_int_sym_potential.rds"))
    fwrite(df_sym_potential_split,  file.path(split_dir, "whole_int_sym_potential.csv"))
    
    # Salva il seed usato per tracciabilità
    cat("Seed usata:", current_seed, "\n", file = file.path(split_dir, "seed_used.txt"))
    
    cat("Seed", current_seed, "| split", n, ": calcolati", length(split_groups),
        "potenziali simmetrici in", split_dir, "\n")
  }
}

################################################################################
### SECONDO BLOCCO: SATURAZIONE DEL POTENZIALE SIMMETRICO WHOLE-INTERFACE
### (da incollare in coda allo script di analisi)
################################################################################

# Cartella prodotta dallo script di split + potenziali
convergence_dir <- "/Users/lorenzosisti/TiNDER/data/whole_int_sym_pot_convergence"
results_dir <- "/Users/lorenzosisti/TiNDER/data/whole_int_sym_pot_convergence_analysis"
dir.create(results_dir, showWarnings = FALSE)

# Numero totale di PDB che hanno prodotto contatti (dal csv dei contatti)
df_contacts <- fread(file.path(convergence_dir, "df_contacts.csv"))
n_total <- uniqueN(df_contacts$pdb_path)

# ---- Leggo i potenziali simmetrici di tutti i seed/split/gruppi ----
rds_files <- list.files(convergence_dir,
                        pattern = "^whole_int_sym_potential\\.rds$",
                        recursive = TRUE, full.names = TRUE)
sym_all <- map_dfr(rds_files, readRDS) %>%
  mutate(pair = paste(resid_i, resid_j, sep = "-"))

# ---- Correlazioni di Pearson a coppie tra i gruppi, per ogni seed e split ----
correlazioni_stats <- sym_all %>%
  group_by(seed, split) %>%
  group_modify(~{
    wide <- pivot_wider(.x, id_cols = pair, names_from = group, values_from = potential)
    cor_mat <- cor(wide[, -1], method = "pearson", use = "pairwise.complete.obs")
    cor_vals <- cor_mat[upper.tri(cor_mat)]
    tibble(
      media_correlazioni = mean(cor_vals, na.rm = TRUE),
      sd_correlazioni = sd(cor_vals, na.rm = TRUE)
    )
  }) %>%
  ungroup() %>%
  rename(split_num = split)

# ---- Media sui seed ----
aggregati <- correlazioni_stats %>%
  group_by(split_num) %>%
  summarise(
    avg_mean = mean(media_correlazioni, na.rm = TRUE),
    avg_sd = mean(sd_correlazioni, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(n_pdb = floor(n_total / split_num))

# Pendenza dei segmenti della curva di saturazione (n_pdb crescente)
aggregati <- aggregati %>%
  arrange(n_pdb) %>%
  mutate(slope = c(NA, diff(avg_mean) / diff(n_pdb)))

print(aggregati)
write.csv(aggregati, file = file.path(results_dir, "aggregati_sym.csv"), row.names = FALSE)

# ---- Tema comune (sfondo bianco) ----
tema_bianco <- theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))

# ---- Grafico di saturazione: media delle correlazioni +- sd media (barre verticali) ----
plot_saturation <- ggplot(aggregati, aes(x = n_pdb, y = avg_mean)) +
  geom_line(color = "#CC4678FF", linewidth = 1) +
  geom_point(color = "#CC4678FF", size = 2) +
  geom_errorbar(aes(ymin = avg_mean - avg_sd, ymax = avg_mean + avg_sd),
                width = diff(range(aggregati$n_pdb)) * 0.01,
                color = "#CC4678FF", alpha = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(
    x = TeX("Number of PDB files used to compute the statistical potentials"),
    y = TeX("Average Pearson's $\\textit{r}$ between potential matrices ($\\pm \\sigma$)")
  ) +
  tema_bianco

print(plot_saturation)

ggsave(file.path(results_dir, "saturation_sym.png"),
       plot = plot_saturation, width = 7, height = 5, dpi = 300, bg = "white")

# ---- plot5: deviazione standard media in funzione del numero di PDB ----
plot5 <- ggplot(aggregati, aes(x = n_pdb, y = avg_sd)) +
  geom_line(color = "#CC4678FF", linewidth = 1) +
  geom_point(color = "#CC4678FF", size = 2) +
  labs(
    x = TeX("Number of PDB files used to compute the statistical potentials"),
    y = TeX("Average $\\sigma$ between potential matrices")
  ) +
  tema_bianco

print(plot5)

ggsave(file.path(results_dir, "avg_sd_sym.png"),
       plot = plot5, width = 7, height = 5, dpi = 300, bg = "white")

# ---- plot6: pendenza dei segmenti del grafico di saturazione ----
plot6 <- ggplot(filter(aggregati, !is.na(slope)), aes(x = n_pdb, y = slope)) +
  geom_line(color = "#CC4678FF", linewidth = 1) +
  geom_point(color = "#CC4678FF", size = 2) +
  labs(
    x = TeX("Number of PDB files used to compute the statistical potentials"),
    y = TeX("Slope of the segments in the convergence plot")
  ) +
  tema_bianco

print(plot6)

ggsave(file.path(results_dir, "slope_sym.png"),
       plot = plot6, width = 7, height = 5, dpi = 300, bg = "white")

ggsave(file.path(results_dir, "saturation_sym.pdf"),
       plot = plot_saturation, width = 8, height = 6, device = cairo_pdf)

ggsave(file.path(results_dir, "avg_sd_sym.pdf"),
       plot = plot5, width = 8, height = 6, device = cairo_pdf)

ggsave(file.path(results_dir, "slope_sym.pdf"),
       plot = plot6, width = 8, height = 6, device = cairo_pdf)

