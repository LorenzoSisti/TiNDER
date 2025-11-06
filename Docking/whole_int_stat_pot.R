# Importing libraries
library(bio3d)
library(dplyr)
library(tidyr)
library(future)
library(furrr)
library(progressr)

### Imposta parallelizzazione per rendere più veloce l'esecuzione del codice
plan(multisession, workers = parallel::detectCores() - 1)
#handlers(global = TRUE)
#handlers("rstudio")

# Setting paths
pdb_dir <- "/Users/lorenzosisti/Downloads/models/"
results_dir <- "/Users/lorenzosisti/Downloads/potenziali_statistici_HDOCK"
dir.create(results_dir, showWarnings = FALSE)

DistCutoff <- 8.5
Nsteps <- 3

pdbL <- list.files(pdb_dir)

# Amino acids ordered according to the Kyte-Doolittle hydrophobicity scale
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE")

### MAIN ###

asym_potentials_whole_int <- read.csv("/Users/lorenzosisti/Downloads/matrice_singola_novembre/V_asym_df.csv")
sym_potentials_whole_int <- read.csv("/Users/lorenzosisti/Downloads/matrice_singola_novembre/V_sym_df.csv")

compute_sym_pot <- function(s) {
  tryCatch({
    
    # Retrieve and read the PDB file
    path_aus <- paste0(pdb_dir, pdbL[s])
    pdb_aus <- read.pdb(path_aus)  
    file_name <- pdbL[s]
    
    parts <- strsplit(file_name, "_")[[1]]
    chain_HL <- c(parts[2], parts[3])
    chain_AG <- parts[4]
    
    # Definiamo un DB solo delle catene laterali (CA per Gly) per poter poi calcolare il relativo centroide
    df_coord <- pdb_aus$atom[pdb_aus$atom$type != "HETATM" & (
      (pdb_aus$atom$resid != "GLY" & !(pdb_aus$atom$elety %in% c("N", "CA", "C", "O"))) |
        (pdb_aus$atom$resid == "GLY" & pdb_aus$atom$elety == "CA")),]
    
    if (nrow(df_coord) == 0) {
      writeLines(paste("File con df_coord vuoto:", path_aus), "errors.log", append = TRUE)
      next
    }
    
    # Eventualmente, rinumeriamo per evitare situazioni del tipo SER100, ASP100A, GLN101. Partiamo con il creare una copia di df_corrected
    df_coord_corrected <- df_coord
    corrected_dfs <- list()
    for (chain in unique(df_coord_corrected$chain)) {
      residui_unici <- unique(paste(df_coord_corrected$resno[df_coord_corrected$chain == chain], 
                                    df_coord_corrected$insert[df_coord_corrected$chain == chain], sep=""))
      nuova_numerazione <- setNames(seq_along(residui_unici), residui_unici)
      df_coord_corrected$resno[df_coord_corrected$chain == chain] <- nuova_numerazione[paste(df_coord_corrected$resno[df_coord_corrected$chain == chain], 
                                                                                             df_coord_corrected$insert[df_coord_corrected$chain == chain], sep="")]
      corrected_dfs[[chain]] <- df_coord_corrected[df_coord_corrected$chain == chain, ]
    }
    
    # Rimettiamo tutto insieme in un unico df
    df_coord_corrected <- do.call(rbind, corrected_dfs)
    
    #Calcoliamo i centroidi per ogni catena laterale
    centroidi_df <- as.data.frame(
      df_coord_corrected %>%
        group_by(chain, resno, resid) %>%
        summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
    )
    
    res_names <- paste(centroidi_df$resid, centroidi_df$resno, centroidi_df$chain, sep = "_")
    df_coord_resid_xyz <- centroidi_df[, c("resid", "resno", "x", "y", "z")] 
    rownames(df_coord_resid_xyz) <- res_names
    
    #Definiamo una matrice delle distanze del DF dei centroidi
    DistMat <- as.matrix(dist(df_coord_resid_xyz[, c("x", "y", "z")]))
    DistMat_bin <- ifelse(DistMat <= DistCutoff, 1, 0) # Matrice delle distanze binaria
    
    condA <- !(centroidi_df$chain %in% chain_HL)
    condHL <- centroidi_df$chain %in% chain_HL
    
    # Definiamo una matrice delle distanze intermolecolari tra i centrodi di Ab e Ag
    Inter_DistMat <- DistMat[condA, condHL]
    Inter_DistMat_Bin <- ifelse(Inter_DistMat <= DistCutoff, 1, 0)
    
    BS_A <- apply(Inter_DistMat_Bin, 1, sum)
    BS_A <- BS_A[BS_A != 0]
    BS_HL <- apply(Inter_DistMat_Bin, 2, sum)
    BS_HL <- BS_HL[BS_HL != 0]
    
    # Identifichiamo i residui facenti parte dei due BS
    BS_A_names <- names(BS_A)
    BS_HL_names <- names(BS_HL)
    
    contact_indices <- which(Inter_DistMat_Bin == 1, arr.ind = TRUE)
    contacts <- data.frame(
      res1 = rownames(Inter_DistMat_Bin)[contact_indices[, 1]],
      res2 = colnames(Inter_DistMat_Bin)[contact_indices[, 2]]
    )
    
    # Estraggo solo l’amminoacido da res1 e res2 (prima del "_")
    #contatti_clean <- contacts %>%
    #  mutate(
    #    aa1 = sub("_.*", "", res1),
    #    aa2 = sub("_.*", "", res2)
    #  ) %>%
    #  count(aa1, aa2, name = "n_contacts")  # Conta quante volte compare la coppia
    
    # Se non ci sono contatti, ritorno NULL
    if (nrow(contacts) == 0) {
      message("⚠️  Nessun contatto trovato in ", pdbL[s])
      return(NULL)
    }
    
    contatti_clean <- contacts %>%
      mutate(
        aa1 = sub("_.*", "", res1),
        aa2 = sub("_.*", "", res2)
      ) %>%
      group_by(aa1, aa2) %>%
      summarise(n_contacts = n(), .groups = "drop")
    
    
    # Merge con potenziali simmetrici
    result <- sym_potentials_whole_int %>%
      select(aa1, aa2, value) %>%
      left_join(contatti_clean, by = c("aa1", "aa2")) %>%
      mutate(
        n_contacts = replace_na(n_contacts, 0),
        potenziale_totale = value * n_contacts
      )
    
    return(result)
    
  })
}

compute_asym_pot <- function(s) {
  tryCatch({
    
    # Retrieve and read the PDB file
    path_aus <- paste0(pdb_dir, pdbL[s])
    pdb_aus <- read.pdb(path_aus)  
    file_name <- pdbL[s]
    
    parts <- strsplit(file_name, "_")[[1]]
    chain_HL <- c(parts[2], parts[3])
    chain_AG <- parts[4]
    
    # Definiamo un DB solo delle catene laterali (CA per Gly) per poter poi calcolare il relativo centroide
    df_coord <- pdb_aus$atom[pdb_aus$atom$type != "HETATM" & (
      (pdb_aus$atom$resid != "GLY" & !(pdb_aus$atom$elety %in% c("N", "CA", "C", "O"))) |
        (pdb_aus$atom$resid == "GLY" & pdb_aus$atom$elety == "CA")),]
    
    if (nrow(df_coord) == 0) {
      writeLines(paste("File con df_coord vuoto:", path_aus), "errors.log", append = TRUE)
      next
    }
    
    # Eventualmente, rinumeriamo per evitare situazioni del tipo SER100, ASP100A, GLN101. Partiamo con il creare una copia di df_corrected
    df_coord_corrected <- df_coord
    corrected_dfs <- list()
    for (chain in unique(df_coord_corrected$chain)) {
      residui_unici <- unique(paste(df_coord_corrected$resno[df_coord_corrected$chain == chain], 
                                    df_coord_corrected$insert[df_coord_corrected$chain == chain], sep=""))
      nuova_numerazione <- setNames(seq_along(residui_unici), residui_unici)
      df_coord_corrected$resno[df_coord_corrected$chain == chain] <- nuova_numerazione[paste(df_coord_corrected$resno[df_coord_corrected$chain == chain], 
                                                                                             df_coord_corrected$insert[df_coord_corrected$chain == chain], sep="")]
      corrected_dfs[[chain]] <- df_coord_corrected[df_coord_corrected$chain == chain, ]
    }
    
    # Rimettiamo tutto insieme in un unico df
    df_coord_corrected <- do.call(rbind, corrected_dfs)
    
    #Calcoliamo i centroidi per ogni catena laterale
    centroidi_df <- as.data.frame(
      df_coord_corrected %>%
        group_by(chain, resno, resid) %>%
        summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
    )
    
    res_names <- paste(centroidi_df$resid, centroidi_df$resno, centroidi_df$chain, sep = "_")
    df_coord_resid_xyz <- centroidi_df[, c("resid", "resno", "x", "y", "z")] 
    rownames(df_coord_resid_xyz) <- res_names
    
    #Definiamo una matrice delle distanze del DF dei centroidi
    DistMat <- as.matrix(dist(df_coord_resid_xyz[, c("x", "y", "z")]))
    DistMat_bin <- ifelse(DistMat <= DistCutoff, 1, 0) # Matrice delle distanze binaria
    
    condA <- !(centroidi_df$chain %in% chain_HL)
    condHL <- centroidi_df$chain %in% chain_HL
    
    # Definiamo una matrice delle distanze intermolecolari tra i centrodi di Ab e Ag
    Inter_DistMat <- DistMat[condA, condHL]
    Inter_DistMat_Bin <- ifelse(Inter_DistMat <= DistCutoff, 1, 0)
    
    BS_A <- apply(Inter_DistMat_Bin, 1, sum)
    BS_A <- BS_A[BS_A != 0]
    BS_HL <- apply(Inter_DistMat_Bin, 2, sum)
    BS_HL <- BS_HL[BS_HL != 0]
    
    # Identifichiamo i residui facenti parte dei due BS
    BS_A_names <- names(BS_A)
    BS_HL_names <- names(BS_HL)
    
    contact_indices <- which(Inter_DistMat_Bin == 1, arr.ind = TRUE)
    contacts <- data.frame(
      res1 = rownames(Inter_DistMat_Bin)[contact_indices[, 1]],
      res2 = colnames(Inter_DistMat_Bin)[contact_indices[, 2]]
    )
    
    #contatti_clean <- contacts %>%
    #  mutate(
    #    aa1 = paste0(sub("_.*", "", res2), "_", ifelse(grepl("_H$|_L$", res2), "Ab", "Ag")),
    #    aa2 = paste0(sub("_.*", "", res1), "_", ifelse(grepl("_H$|_L$", res1), "Ab", "Ag"))
   #   ) %>%
   #   count(aa1, aa2, name = "n_contacts")
    
    # Se non ci sono contatti, ritorno NULL
    if (nrow(contacts) == 0) {
      message("⚠️  Nessun contatto trovato in ", pdbL[s])
      return(NULL)
    }
    
    contatti_clean <- contacts %>%
      mutate(
        aa1 = paste0(sub("_.*", "", res2), "_", ifelse(grepl("_H$|_L$", res2), "Ab", "Ag")),
        aa2 = paste0(sub("_.*", "", res1), "_", ifelse(grepl("_H$|_L$", res1), "Ab", "Ag"))
      ) %>%
      group_by(aa1, aa2) %>%
      summarise(n_contacts = n(), .groups = "drop")
    
    
    result <- asym_potentials_whole_int %>%
      select(aa1, aa2, value) %>%
      left_join(contatti_clean, by = c("aa1", "aa2")) %>%
      mutate(
        n_contacts = replace_na(n_contacts, 0),
        potenziale_totale = value * n_contacts
      )
    
    return(result)
    
  })
}

### RUN SU TUTTI I FILE PDB ###

summary_results <- future_map_dfr(
  seq_along(pdbL),
  function(i) {
    message("Processing: ", pdbL[i])
    
    res_sym  <- compute_sym_pot(i)
    res_asym <- compute_asym_pot(i)
    
    if (is.null(res_sym) | is.null(res_asym)) {
      return(data.frame(
        pdb = pdbL[i],
        sum_sym = NA,
        sum_asym = NA,
        mean_sym = NA,
        mean_asym = NA
      ))
    }
    
    data.frame(
      pdb = pdbL[i],
      sum_sym  = sum(res_sym$potenziale_totale, na.rm = TRUE),
      sum_asym = sum(res_asym$potenziale_totale, na.rm = TRUE),
      mean_sym  = mean(res_sym$potenziale_totale[res_sym$n_contacts > 0], na.rm = TRUE),
      mean_asym = mean(res_asym$potenziale_totale[res_asym$n_contacts > 0], na.rm = TRUE)
    )
  },
  .progress = TRUE   
)


# Salvataggio del CSV
write.csv(summary_results,
          file = file.path(results_dir, "potenziali_per_modello.csv"),
          row.names = FALSE)

message("\n✅ Analisi completata!")
message("➡️ Risultato salvato in: ", file.path(results_dir, "potenziali_per_modello.csv"))
