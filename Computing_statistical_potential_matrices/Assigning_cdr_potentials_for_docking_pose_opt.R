# Caricamento delle librerie necessarie (aggiunto readr per read_csv)
pacman::p_load(bio3d, 
               reshape2,
               ftrCOOL,
               AATtools,
               ggplot2,
               readr,
               tidyr,
               purrr,
               dplyr)

# Cartella contenente i file PDB
pdb_dir <- "/Users/lorenzosisti/Downloads/docked_structures_renamed_AF3_11_06"
results_dir <- "/Users/lorenzosisti/Downloads/dock_potenziali_statistici_cdr_12_06/"
dir.create(results_dir, showWarnings = FALSE)

pdb_files <- list.files(path = pdb_dir, pattern = "\\.pdb$", full.names = TRUE)

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

# Caricamento dei potenziali
sym_cdr_potentials  <- read_csv("/Users/lorenzosisti/Downloads/potenziali_statistici_cdr_12_06/V_sym_wide.csv", show_col_types = FALSE)
asym_cdr_potentials <- read_csv("/Users/lorenzosisti/Downloads/potenziali_statistici_cdr_12_06/V_asym_wide.csv", show_col_types = FALSE)

# Trasformazione ASIMMETRICI in formato Long
asym_long <- asym_cdr_potentials %>%
  pivot_longer(
    cols = starts_with("V_"), 
    names_to = "cdr", 
    names_prefix = "V_", 
    values_to = "potential_asym"
  )

# Trasformazione SIMMETRICI in formato Long
sym_long <- sym_cdr_potentials %>%
  pivot_longer(
    cols = starts_with("V_"), 
    names_to = "cdr", 
    names_prefix = "V_", 
    values_to = "potential_sym"
  )

# Inizializzazione della lista master per i risultati finali per PDB
pdb_scores_list <- list()

#pdb_file <- "/Users/lorenzosisti/Downloads/database_settembre_renamed/6wzm.pdb" 

# Inizio Loop su tutti i file PDB  
for (pdb_file in pdb_files) { 
  
  tryCatch({
    
    # Lettura della struttura PDB  
    pdb_aus <- read.pdb(pdb_file) 
    chain_HL <- c("H", "L")
    chain_AG <- "A" 
    
    # Coarse-grained representation
    df_coord <- pdb_aus$atom[pdb_aus$atom$type != "HETATM" & (
      (pdb_aus$atom$resid != "GLY" & !(pdb_aus$atom$elety %in% c("N", "CA", "C", "O"))) |
        (pdb_aus$atom$resid == "GLY" & pdb_aus$atom$elety == "CA")),]
    
    if (nrow(df_coord) == 0) {
      cat("File con df_coord vuoto:", pdb_file, "\n", file = "errors.log", append = TRUE)
      next 
    }
    
    df_coord_corrected <- df_coord
    
    # 1. Gestione dei valori NA nella colonna insert
    df_coord_corrected$insert[is.na(df_coord_corrected$insert)] <- ""
    
    # 2. Compute centroid coordinates
    centroidi_df <- as.data.frame(
      df_coord_corrected %>%
        group_by(chain, resno, insert, resid) %>%
        summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
    )
    
    # 3. Creazione di nomi univoci
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
    condA <- !(centroidi_df$chain %in% c("H", "L"))
    condHL <- centroidi_df$chain %in% c("H", "L")   
    
    # Subset distance matrix (Righe = Antigene, Colonne = Anticorpo)
    Inter_DistMat_Bin <- DistMat_bin[condA, condHL, drop = FALSE]
    
    # Vettore per tracciare tutti i residui dell'antigene che toccano una qualsiasi CDR
    ag_interacting_with_any_cdr <- c()
    
    # ---------------------------------------------------------
    # CORREZIONE 1: Inizializzazione lista contatti per QUESTO PDB
    # ---------------------------------------------------------
    current_pdb_contacts <- list() 
    contact_counter <- 1
    
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
        next 
      }
      
      Inter_DistMat_Bin_CDR <- Inter_DistMat_Bin[, valid_cols, drop = FALSE]
      
      # Calcola i contatti specifici per questa CDR
      BS_A_cdr <- apply(Inter_DistMat_Bin_CDR, 1, sum)
      BS_A_cdr <- BS_A_cdr[BS_A_cdr != 0]
      
      BS_HL_cdr <- apply(Inter_DistMat_Bin_CDR, 2, sum)
      BS_HL_cdr <- BS_HL_cdr[BS_HL_cdr != 0]
      
      BS_A_names_cdr <- names(BS_A_cdr)
      BS_HL_names_cdr <- names(BS_HL_cdr)
      
      # Aggiunge i residui Ag all'ensemble globale di questo PDB
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
            
            # Matrice Asimmetrica Cumulative
            contact_matrices_asym_rings_sum[[cdr_name]][aa_ab, aa_ag] <- contact_matrices_asym_rings_sum[[cdr_name]][aa_ab, aa_ag] + 1
            
            # Matrice Simmetrica Cumulative
            if (aa_ab == aa_ag) {
              contact_matrices_sym_rings_sum[[cdr_name]][aa_ab, aa_ag] <- contact_matrices_sym_rings_sum[[cdr_name]][aa_ab, aa_ag] + 1
            } else {
              contact_matrices_sym_rings_sum[[cdr_name]][aa_ab, aa_ag] <- contact_matrices_sym_rings_sum[[cdr_name]][aa_ab, aa_ag] + 1
              contact_matrices_sym_rings_sum[[cdr_name]][aa_ag, aa_ab] <- contact_matrices_sym_rings_sum[[cdr_name]][aa_ag, aa_ab] + 1
            }
            
            # ---------------------------------------------------------
            # CORREZIONE 2: Salvataggio effettivo del singolo contatto
            # ---------------------------------------------------------
            current_pdb_contacts[[contact_counter]] <- data.frame(
              cdr = cdr_name,
              aa1 = aa_ab,  # Anticorpo
              aa2 = aa_ag,  # Antigene
              stringsAsFactors = FALSE
            )
            contact_counter <- contact_counter + 1
          }
        }
      }
    } # Fine iterazione CDR
    
    # Aggiornamento residui antigene cumulativi (ripristinato dal codice originale)
    if (length(ag_interacting_with_any_cdr) > 0) {
      df_centroidi_ag_bs <- centroidi_df[rownames(centroidi_df) %in% ag_interacting_with_any_cdr, ]
      residue_counts_interface_ligand_sum <- residue_counts_interface_ligand_sum + table(factor(df_centroidi_ag_bs$resid, levels = amino_acids))
    }
    
    # =========================================================
    # COLLASSO E CALCOLO POTENZIALI (ASYM e SYM) PER IL PDB
    # =========================================================
    if (length(current_pdb_contacts) > 0) {
      
      # 1. Unisce tutti i singoli contatti tracciati nel PDB
      df_contacts <- bind_rows(current_pdb_contacts)
      
      # 2. Raggruppa e conta le frequenze dei contatti
      df_scored <- df_contacts %>%
        group_by(cdr, aa1, aa2) %>%
        summarise(count = n(), .groups = "drop") %>%
        
        # 3. Join con i potenziali Asimmetrici
        left_join(asym_long, by = c("cdr", "aa1", "aa2")) %>%
        
        # 4. Join con i potenziali Simmetrici
        left_join(sym_long, by = c("cdr", "aa1", "aa2")) %>%
        
        # 5. Calcolo dei potenziali totali (moltiplicati per il numero di contatti)
        mutate(
          total_potential_asym = count * potential_asym,
          total_potential_sym  = count * potential_sym
        )
      
      # 6. Salvataggio nella lista master usando il nome del file PDB
      pdb_name <- basename(pdb_file)
      pdb_scores_list[[pdb_name]] <- df_scored
      
    }
    
  }, error = function(e) {
    cat(paste("Errore elaborazione file:", pdb_file, "-", e$message, "\n"), file = "errors.log", append = TRUE)
  })
}

# 1. Collassa la lista di tutti i PDB in un unico grande Dataframe
final_global_dataframe <- bind_rows(pdb_scores_list, .id = "pdb_filename")

# 2. Salva il file nella cartella dei risultati
write_csv(final_global_dataframe, 
          file = file.path(results_dir, "tutti_i_potenziali_singoli_PDB.csv"))

cat("Salvataggio completato! Trovi il file master in:", results_dir, "\n")