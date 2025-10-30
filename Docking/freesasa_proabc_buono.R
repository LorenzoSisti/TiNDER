# ------- CONFIGURAZIONE -------
base_dir <- "/Users/lorenzosisti/Downloads/proabc2-prediction"
output_dir <- file.path(base_dir, "combined_chothia_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# mapping 1-letter -> 3-letter
aa1to3 <- c(
  A = "ALA", R = "ARG", N = "ASN", D = "ASP", C = "CYS",
  Q = "GLN", E = "GLU", G = "GLY", H = "HIS", I = "ILE",
  L = "LEU", K = "LYS", M = "MET", F = "PHE", P = "PRO",
  S = "SER", T = "THR", W = "TRP", Y = "TYR", V = "VAL"
)

# ------- trova file heavy e light -------
heavy_files <- list.files(path = base_dir, pattern = "heavy-pred\\.csv$", recursive = TRUE, full.names = TRUE)
light_files <- list.files(path = base_dir, pattern = "light-pred\\.csv$", recursive = TRUE, full.names = TRUE)

# ottieni cartelle uniche
folders <- unique(dirname(c(heavy_files, light_files)))

# ------- funzione compatta per processare dataframe -------
process_and_tag <- function(df, chain_label) {
  df_sel <- df[df$pt >= 0.4, , drop = FALSE]
  if (nrow(df_sel) == 0) return(df_sel)
  df_sel$Sequence <- toupper(df_sel$Sequence)
  df_sel$ResName3 <- aa1to3[df_sel$Sequence]
  df_sel$ResName3[is.na(df_sel$ResName3)] <- "UNK"
  df_sel$Residue <- paste(df_sel$ResName3, df_sel$Chothia, sep = "_")
  df_sel$Chain <- chain_label
  return(df_sel)
}

# ------- loop principale e raccolta globale -------
all_combined <- list()   # lista per il CSV globale

for (folder in folders) {
  pdb_id <- basename(folder)
  combined <- list()
  
  heavy_path <- file.path(folder, "heavy-pred.csv")
  if (file.exists(heavy_path)) {
    df_h <- read.csv(heavy_path, stringsAsFactors = FALSE)
    processed_h <- process_and_tag(df_h, "H")
    if (nrow(processed_h) > 0) combined[[length(combined) + 1]] <- processed_h
  }
  
  light_path <- file.path(folder, "light-pred.csv")
  if (file.exists(light_path)) {
    df_l <- read.csv(light_path, stringsAsFactors = FALSE)
    processed_l <- process_and_tag(df_l, "L")
    if (nrow(processed_l) > 0) combined[[length(combined) + 1]] <- processed_l
  }
  
  if (length(combined) > 0) {
    combined_df <- do.call(rbind, combined)
    # aggiungi colonna PDB per il file globale
    combined_df$PDB <- pdb_id
    
    # riordina colonne (opzionale)
    #front <- c("PDB", "Chain", "ResName3", "Chothia", "Residue", "pt")
    #present_front <- intersect(front, names(combined_df))
    #rest <- setdiff(names(combined_df), present_front)
    #combined_df <- combined_df[c(present_front, rest)]
    
    # salva file singolo per PDB
    out_file <- file.path(output_dir, paste0(pdb_id, "_chothia.csv"))
    write.csv(combined_df, out_file, row.names = FALSE)
    message("Saved: ", out_file, " (rows: ", nrow(combined_df), ")")
    
    # aggiungi al raccolto globale
    all_combined[[length(all_combined) + 1]] <- combined_df
  }
}

# ------- salva il CSV globale se ci sono dati -------
if (length(all_combined) > 0) {
  global_df <- do.call(rbind, all_combined)
  global_out <- file.path(output_dir, "all_chothia_combined.csv")
  write.csv(global_df, global_out, row.names = FALSE)
  message("Saved global combined CSV: ", global_out, " (rows: ", nrow(global_df), ")")
} else {
  message("No records found to combine globally.")
}

global_proabc_df <- read.csv("/Users/lorenzosisti/Downloads/proabc2-prediction/combined_chothia_results/all_chothia_combined.csv")
global_proabc_df$ResID <- toupper(paste(global_proabc_df$PDB, global_proabc_df$Residue, global_proabc_df$Chain, sep = "_"))
#Ora creo un ponte per fare un merge con freeSASA

global_proabc_df[global_proabc_df$PDB == "8trt", ]

freesasa_dir <- "/Users/lorenzosisti/Downloads/SASA_for_docking"
sasa_file <- read.csv(paste(freesasa_dir, "/combined_sasa_rsa.csv", sep = ""))
sasa_file$ResID <- paste0(
  trimws(toupper(sub("_chothia\\.pdb$", "", as.character(sasa_file$PDB), ignore.case = TRUE))), "_",
  trimws(as.character(sasa_file$ResName)), "_",
  trimws(as.character(sasa_file$ResNum)),  "_",
  trimws(as.character(sasa_file$Chain))
)

# --- 1) Normalizza PDB e colonne (Il tuo codice qui è corretto) ---
#sasa_file$PDB_clean <- sub("_chothia\\.pdb$", "", as.character(sasa_file$PDB), ignore.case = TRUE)
#sasa_file$PDB_clean <- toupper(trimws(as.character(sasa_file$PDB_clean)))
#sasa_file$Chain       <- toupper(trimws(as.character(sasa_file$Chain)))
#sasa_file$ResNum      <- trimws(as.character(sasa_file$ResNum)) 


#global_proabc_df$PDB     <- toupper(trimws(as.character(global_proabc_df$PDB)))
#global_proabc_df$Chain   <- toupper(trimws(as.character(global_proabc_df$Chain)))
#global_proabc_df$Chothia <- trimws(as.character(global_proabc_df$Chothia))

# --- 2) Fai il merge AGGIUNGENDO i suffissi per colonne duplicate (es. ResID) ---
merged_data <- merge(
  global_proabc_df, # Metti proabc per primo per mantenere le sue colonne
  sasa_file,
  by.x = "ResID",
  by.y = "ResID",
  all = FALSE,  # inner join -> solo intersezione
  suffixes = c("_proabc", "_sasa") # Gestisce ResID, PDB, ecc. duplicati
)

# --- 3) APPLICA IL FILTRO RSA > 0.25 ---
# Ora merged_data contiene le colonne di entrambi; possiamo filtrare usando la colonna RSA
proabc_filtered_by_rsa <- merged_data[which(merged_data$RSA > 0.25), ]

# controlli rapidi
message("Totale righe ProABC prima: ", nrow(global_proabc_df))
message("Totale matches (sasa ∩ proabc): ", nrow(merged_data))
message("Matches con RSA > 0.25: ", nrow(proabc_filtered_by_rsa))

# --- 4) Salva i dati PROABC filtrati ---
# Selezioniamo le colonne originali di proabc (e magari RSA per controllo)
proabc_original_cols <- names(global_proabc_df)
cols_to_save <- c(proabc_original_cols, "RSA") # Aggiungiamo RSA

# Assicurati che tutte le colonne esistano nel df unito
# (le chiavi di merge by.x "PDB", "Chain", "Chothia" sono già in proabc_original_cols)
final_cols <- intersect(cols_to_save, names(proabc_filtered_by_rsa))

proabc_final_output <- proabc_filtered_by_rsa[, final_cols]

# Salva il dataframe proabc filtrato
out_file <- file.path(output_dir, "proabc_filtered_by_RSA_gt025.csv")
write.csv(proabc_final_output, out_file, row.names = FALSE)

message("Salvato file proabc filtrato per RSA: ", out_file)









# --- 1) Normalizza PDB e colonne come carattere e senza spazi ---
# rimuovi suffix tipo "_chothia.pdb" dai PDB di freesasa
sasa_file$PDB_clean <- sub("_chothia\\.pdb$", "", as.character(sasa_file$PDB), ignore.case = TRUE)

# trim e tocca tutto a maiuscolo per evitare mismatch case/whitespace
sasa_file$PDB_clean <- toupper(trimws(as.character(sasa_file$PDB_clean)))
sasa_file$Chain     <- toupper(trimws(as.character(sasa_file$Chain)))
sasa_file$ResNum    <- trimws(as.character(sasa_file$ResNum))   # conserva "100A" così com'è

global_proabc_df$PDB   <- toupper(trimws(as.character(global_proabc_df$PDB)))
global_proabc_df$Chain <- toupper(trimws(as.character(global_proabc_df$Chain)))
# Chothia può contenere "100A" ecc. -> trattalo come stringa
global_proabc_df$Chothia <- trimws(as.character(global_proabc_df$Chothia))

# --- 2) fai il merge ESATTO su (PDB, Chain, ResNum/Chothia) ---
merged_exact <- merge(
  sasa_file, global_proabc_df,
  by.x = c("PDB_clean", "Chain", "ResNum"),
  by.y = c("PDB",       "Chain", "Chothia"),
  all = FALSE    # inner join -> solo intersezione
)

# se desideri solo le righe di sasa_file filtrate (senza colonne extra di proabc)
sasa_matched <- merged_exact[, intersect(names(sasa_file), names(merged_exact)), drop = FALSE]

# controlli rapidi
message("Totale sasa_file prima: ", nrow(sasa_file))
message("Totale matches trovati (sasa ∩ proabc): ", nrow(sasa_matched))

# salva se vuoi
write.csv(sasa_matched, file.path(output_dir, "sasa_matched_to_proabc_byResNum.csv"), row.names = FALSE)



















# --- normalizzazione e preparazione chiavi ---
# togli suffix come "_chothia.pdb" da sasa_file$PDB per ottenere l'ID PDB compatibile
sasa_file$PDB_clean <- sub("_chothia\\.pdb$", "", as.character(sasa_file$PDB), ignore.case = TRUE)
sasa_file$PDB_clean <- trimws(sasa_file$PDB_clean)
sasa_file$Chain <- trimws(as.character(sasa_file$Chain))
sasa_file$ResNum <- trimws(as.character(sasa_file$ResNum))
# assicurati che il numero di residuo sia numerico (o carattere coerente) nella colonna usata per il merge
#sasa_file$Chothia <- as.integer(sasa_file$ResNum)   # crea colonna 'Chothia' per coerenza con global_proabc_df
sasa_file$Chothia <- sasa_file$ResNum   # crea colonna 'Chothia' per coerenza con global_proabc_df

# normalizza global_proabc_df
global_proabc_df$PDB <- trimws(as.character(global_proabc_df$PDB))
global_proabc_df$Chain <- trimws(as.character(global_proabc_df$Chain))
#global_proabc_df$Chothia <- as.integer(as.character(global_proabc_df$Chothia)) # già presente ma assicuriamolo

# opzionale: rendi PDB tutti minuscoli o maiuscoli per evitare mismatch di case
sasa_file$PDB_clean <- toupper(sasa_file$PDB_clean)
global_proabc_df$PDB <- toupper(global_proabc_df$PDB)

# --- INNER JOIN (solo righe presenti in entrambi) ---
merged <- merge(
  sasa_file, global_proabc_df,
  by.x = c("PDB_clean", "Chain", "Chothia"),
  by.y = c("PDB",       "Chain", "Chothia"),
  all = FALSE    # inner join
)

# merged contiene colonne di entrambe le tabelle; se vuoi solo le righe di freesasa filtrate:
sasa_matched <- merged[, intersect(names(sasa_file), names(merged)), drop = FALSE]

# oppure, se preferisci mantenere tutte le colonne di sasa_file più alcune scelte da global_proabc_df:
# sasa_matched2 <- merged[, c(names(sasa_file), "Sequence", "pt", "ResName3"), drop = FALSE]

# verifica
message("Righe originali sasa_file: ", nrow(sasa_file))
message("Righe global_proabc_df (candidati): ", nrow(global_proabc_df))
message("Righe matched (sasa ∩ proabc): ", nrow(sasa_matched))

# salva risultato se vuoi
write.csv(sasa_matched, file.path(output_dir, "sasa_matched_to_proabc.csv"), row.names = FALSE)

