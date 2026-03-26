library(dplyr)
library(purrr)
library(stringr)

# 1. Cartella GIUSTA con i file CSV del Blocco 6
results_dir <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali_marzo/risultati_scoring_finale_HDOCK_WholeInt/"

# 2. Pattern corretto: legge SOLO i file che iniziano con "scores_" e finiscono in ".csv"
all_csv_files <- list.files(results_dir, pattern = "^scores_.*\\.csv$", full.names = TRUE)

cat("Trovati", length(all_csv_files), "file CSV da unire...\n")

# Leggi tutti i file e combinali
full_results_df <- map_dfr(all_csv_files, function(file_path) {
  
  # Estrai info dal nome del file
  filename <- basename(file_path)
  
  # es: "scores_seed_101_split_2_group_1_sym.csv"
  info <- str_match(filename, "scores_(seed_\\d+)_(split_\\d+)_(group_\\d+)_(sym|asym)\\.csv")
  
  # Leggi il CSV
  read.csv(file_path) %>%
    mutate(
      seed = info[, 2],
      split_type = info[, 3],
      group_num = info[, 4],
      potential_type = info[, 5]
    )
})

# Salva il file master
write.csv(full_results_df, file.path(results_dir, "MASTER_RESULTS_HDOCK_tidy.csv"), row.names = FALSE)
saveRDS(full_results_df, file.path(results_dir, "MASTER_RESULTS_HDOCK_tidy.rds"))

cat("Fatto! Risultati combinati salvati.\n")
