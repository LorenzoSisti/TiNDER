### THIS SCRIPTS FILTERS OUT THE DOCKED COMPLEXES WHERE THE ANTIGEN CONTACTS AMINO ACIDS EXTERNAL TO THE CDRs ###

### Required libraries
pacman::p_load(bio3d, dplyr, future, furrr, purrr, progressr,
               pheatmap, patchwork, ggplotify, reshape2, tidyr, data.table, ggplot2)

# Define the path to a custom function files
source("/Users/lorenzosisti/Documents/Script_ottimizzati_funzioni/functions.R")

# Distance cutoff (Å) to define contact between side-chains centroids
DistCutoff <- 8.5  

### Set up parallelization to speed up computation
plan(multisession, workers = parallel::detectCores() - 1)
handlers(global = TRUE)
handlers("rstudio")

pdb_dir <- "/Users/lorenzosisti/Downloads/docked_structures_renamed_AF3_11_06/"
results_dir <- "/Users/lorenzosisti/Downloads/cdrs_filtered_docked_structure/"
dir.create(results_dir, showWarnings = FALSE)

all_pdbs <- list.files(pdb_dir, pattern = "*.pdb", recursive = TRUE, full.names = TRUE)

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
    dt_centroids <- dt_coord[, .(x = mean(x), y = mean(y), z = mean(z)), by = .(chain, resno, resid, region)]
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
        region_ab = region,
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

failed_results <- discard(results_list, ~ .x$ok)
errors_df <- map_dfr(failed_results, ~ tibble(file = .x$filename, error = .x$error))
print(errors_df, n = 50)

### Combina tutti i contatti
df_contacts <- rbindlist(
  map(valid_results, "contacts"),
  use.names = TRUE
)

### Salvataggio
saveRDS(df_contacts,      file.path(results_dir, "df_contacts.rds"))
fwrite(df_contacts,       file.path(results_dir, "df_contacts.csv"))

