### THIS SCRIPT ANALYZES ANTIGEN–ANTIBODY INTERACTIONS FROM MULTIPLE PDB FILES ###
### IT CALCULATES RESIDUE–RESIDUE CONTACT MATRICES AND AMINO ACID COUNTS WITHIN CONCENTRIC "RINGS" AROUND THE INTERFACE CENTER ###
### RESULTS ARE SAVED AS CSV AND RDS FILES FOR EACH RING ###

### Libraries ###
library(bio3d)
library(dplyr)

### Setting paths, directories and global variables ###
pdb_dir <- "/Users/lorenzosisti/Downloads/database_settembre_renamed/"
results_dir_intra <- "/Users/lorenzosisti/Downloads/matrici_stratificate/" 
dir.create(results_dir_intra, showWarnings = FALSE)

DistCutoff <- 8.5 # Two residues are defined as "in contact" if the centroids of their side chain are at less than 8.5 Angstroms
Nsteps <- 3 # Number of concentric "rings" around the interface

pdbL <- list.files(pdb_dir)

# Amino acids ordered according to the Kyte-Doolittle hydrophobicity scale
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE")

### MAIN: GENERATING CONTACT MATRICES AND COMPUTING RESIDUE COUNTS ###

# MODIFICATO: Inizializzazione delle strutture dati cumulative
# Definiamo i nomi per le 6 combinazioni di ring
ring_pairs <- c("1-1", "2-2", "3-3", "1-2", "2-3", "1-3")

# Inizializza le strutture dati CUMULATIVE per le 6 matrici
contact_matrices_asym_rings_sum <- vector("list", length(ring_pairs))
names(contact_matrices_asym_rings_sum) <- ring_pairs
contact_matrices_sym_rings_sum  <- vector("list", length(ring_pairs))
names(contact_matrices_sym_rings_sum) <- ring_pairs

# I conteggi dei residui restano basati sui singoli ring (Nsteps = 3)
residue_counts_interface_ab_sum <- vector("list", Nsteps) 
residue_counts_interface_ligand_df_sum <- vector("list", Nsteps) 

# Matrice di base
base_matrix <- matrix(0, nrow = 20, ncol = 20, dimnames = list(amino_acids, amino_acids))

# Inizializza le 6 matrici cumulative
for (pair in ring_pairs) {
  contact_matrices_asym_rings_sum[[pair]] <- base_matrix
  contact_matrices_sym_rings_sum[[pair]]  <- base_matrix
}

# Inizializza i conteggi per i 3 ring (questo non cambia)
for (step in 1:Nsteps) {
  residue_counts_interface_ab_sum[[step]] <- setNames(rep(0, length(amino_acids)), amino_acids)
  residue_counts_interface_ligand_df_sum[[step]] <- setNames(rep(0, length(amino_acids)), amino_acids)
}


#Function to process each PDB file in the pdbL list
process_pdb_file <- function(s) {
  tryCatch({
    cat("Elaborating:", pdbL[s], "\n")
    
    # Retrieve and read the PDB file
    path_aus <- paste0(pdb_dir, pdbL[s])
    pdb_aus <- read.pdb(path_aus) 
    file_name <- pdbL[s]
    chain_HL <- c("H", "L")
    chain_AG <- "A"
    
    # To speed up contact calculations, we use a coarse-grained representation:
    # ... (codice invariato) ...
    df_coord <- pdb_aus$atom[pdb_aus$atom$type != "HETATM" & (
      (pdb_aus$atom$resid != "GLY" & !(pdb_aus$atom$elety %in% c("N", "CA", "C", "O"))) |
        (pdb_aus$atom$resid == "GLY" & pdb_aus$atom$elety == "CA")),]
    
    if (nrow(df_coord) == 0) {
      cat("File con df_coord vuoto:", path_aus, "\n", file = "errors.log", append = TRUE)
      return(list(ok = FALSE, filename = pdbL[s], error = "df_coord vuoto"))
    }
    
    # ... (codice invariato per correzione numerazione) ...
    df_coord_corrected <- df_coord
    corrected_dfs <- list()
    
    for (chain in unique(df_coord_corrected$chain)) {
      unique_residues <- unique(paste(df_coord_corrected$resno[df_coord_corrected$chain == chain],
                                      df_coord_corrected$insert[df_coord_corrected$chain == chain], sep=""))
      new_numbering <- setNames(seq_along(unique_residues), unique_residues)
      df_coord_corrected$resno[df_coord_corrected$chain == chain] <- new_numbering[paste(df_coord_corrected$resno[df_coord_corrected$chain == chain],
                                                                                         df_coord_corrected$insert[df_coord_corrected$chain == chain], sep="")]
      corrected_dfs[[chain]] <- df_coord_corrected[df_coord_corrected$chain == chain, ]
    }
    
    df_coord_corrected <- do.call(rbind, corrected_dfs)
    
    # ... (codice invariato per calcolo centroidi) ...
    centroids_df <- as.data.frame(
      df_coord_corrected %>%
        group_by(chain, resno, resid) %>%
        summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
    )
    
    res_names <- paste(centroids_df$resid, centroids_df$resno, centroids_df$chain, sep = "_")
    df_coord_resid_xyz <- centroids_df[, c("resid", "resno", "x", "y", "z")]
    rownames(df_coord_resid_xyz) <- res_names
    
    DistMat <- as.matrix(dist(df_coord_resid_xyz[, c("x", "y", "z")])) 
    
    # ... (codice invariato per definizione interfaccia e ring) ...
    condA <- !(centroids_df$chain %in% chain_HL)
    condHL <- centroids_df$chain %in% chain_HL
    
    Inter_DistMat <- DistMat[condA, condHL]
    Inter_DistMat_Bin <- ifelse(Inter_DistMat <= DistCutoff, 1, 0)
    
    BS_A <- apply(Inter_DistMat_Bin, 1, sum)
    BS_A <- BS_A[BS_A != 0]
    BS_HL <- apply(Inter_DistMat_Bin, 2, sum)
    BS_HL <- BS_HL[BS_HL != 0]
    
    if (length(BS_A) == 0 || length(BS_HL) == 0) return(NULL)
    
    BS_A_names <- names(BS_A)
    BS_HL_names <- names(BS_HL)
    
    df_centroids_BS <- df_coord_resid_xyz[rownames(df_coord_resid_xyz) %in% c(BS_A_names, BS_HL_names), ]
    center_BS <- colMeans(df_centroids_BS[, c("x", "y", "z")])
    DistMat_centroid <- as.matrix(dist(rbind(df_centroids_BS[, c("x", "y", "z")], center_BS)))
    vet_dist_centroid <- DistMat_centroid[-nrow(DistMat_centroid), ncol(DistMat_centroid)]
    r_max <- max(vet_dist_centroid)
    borders <- c(0, 0.37, 0.64, 1) * r_max 
    
    # ... (codice invariato per assegnazione contatti e ring) ...
    contact_indices <- which(Inter_DistMat_Bin == 1, arr.ind = TRUE)
    contacts <- data.frame(
      res1 = rownames(Inter_DistMat_Bin)[contact_indices[, 1]],
      res2 = colnames(Inter_DistMat_Bin)[contact_indices[, 2]]
    )
    
    contacts$dist_res1_centroid <- vet_dist_centroid[contacts$res1]
    contacts$dist_res2_centroid <- vet_dist_centroid[contacts$res2]
    
    assign_ring <- function(dist) findInterval(dist, borders, left.open = TRUE)
    contacts$ring_res1 <- sapply(contacts$dist_res1_centroid, assign_ring)
    contacts$ring_res2 <- sapply(contacts$dist_res2_centroid, assign_ring)
    
    
    # MODIFICATO: Inizializzazione delle strutture dati locali
    # Usiamo gli stessi nomi definiti globalmente
    ring_pairs <- c("1-1", "2-2", "3-3", "1-2", "2-3", "1-3")
    
    contact_matrices_asym_rings <- vector("list", length(ring_pairs))
    names(contact_matrices_asym_rings) <- ring_pairs
    contact_matrices_sym_rings  <- vector("list", length(ring_pairs))
    names(contact_matrices_sym_rings) <- ring_pairs
    
    # Matrice di base locale
    local_base_matrix <- matrix(0, nrow=length(amino_acids), ncol=length(amino_acids), dimnames=list(amino_acids, amino_acids))
    
    for (pair in ring_pairs) {
      contact_matrices_asym_rings[[pair]] <- local_base_matrix
      contact_matrices_sym_rings[[pair]]  <- local_base_matrix
    }
    
    # L'inizializzazione dei conteggi resta invariata (per 3 ring)
    local_ab_counts <- vector("list", Nsteps)
    local_ag_counts <- vector("list", Nsteps)
    for (step in 1:Nsteps) {
      local_ab_counts[[step]] <- setNames(rep(0, length(amino_acids)), amino_acids)
      local_ag_counts[[step]] <- setNames(rep(0, length(amino_acids)), amino_acids)
    }
    
    
    # MODIFICATO: Logica di riempimento delle matrici di contatto
    for (k in seq_len(nrow(contacts))) {
      aa1 <- strsplit(contacts$res1[k], "_")[[1]][1]
      aa2 <- strsplit(contacts$res2[k], "_")[[1]][1]
      if (!(aa1 %in% amino_acids && aa2 %in% amino_acids)) next
      
      ring1 <- contacts$ring_res1[k]
      ring2 <- contacts$ring_res2[k]
      
      # NUOVO: Controlla se i ring sono validi (1, 2, o 3)
      if (ring1 %in% 1:Nsteps && ring2 %in% 1:Nsteps) {
        
        # NUOVO: Crea una chiave standardizzata per la coppia di ring
        # Questo ordina i ring, così (1, 2) e (2, 1) diventano entrambi "1-2"
        rings_sorted <- sort(c(ring1, ring2))
        pair_key <- paste(rings_sorted[1], rings_sorted[2], sep = "-")
        
        # La logica sottostante per determinare le catene rimane la stessa
        chain1 <- strsplit(contacts$res1[k], "_")[[1]][3]
        chain2 <- strsplit(contacts$res2[k], "_")[[1]][3]
        
        ## ---- ASYMMETRIC CONTACT MATRICES ----
        # MODIFICATO: Usa 'pair_key' invece di 'assigned_ring'
        
        if (chain1 %in% chain_HL && chain2 %in% chain_AG) {
          # antibody → antigen
          contact_matrices_asym_rings[[pair_key]][aa1, aa2] <- contact_matrices_asym_rings[[pair_key]][aa1, aa2] + 2
        } else if (chain1 %in% chain_AG && chain2 %in% chain_HL) {
          # antigen → antibody (switch to keep rows = Ab, columns = Ag)
          contact_matrices_asym_rings[[pair_key]][aa2, aa1] <- contact_matrices_asym_rings[[pair_key]][aa2, aa1] + 2
        }
        
        ## ---- SYMMETRIC CONTACT MATRICES ----
        # MODIFICATO: Usa 'pair_key' invece di 'assigned_ring'
        
        if (aa1 == aa2) {
          contact_matrices_sym_rings[[pair_key]][aa1, aa1] <- contact_matrices_sym_rings[[pair_key]][aa1, aa1] + 2
        } else {
          contact_matrices_sym_rings[[pair_key]][aa1, aa2] <- contact_matrices_sym_rings[[pair_key]][aa1, aa2] + 2
          contact_matrices_sym_rings[[pair_key]][aa2, aa1] <- contact_matrices_sym_rings[[pair_key]][aa2, aa1] + 2
        }
      } # Fine del blocco if (ring1 == ring2 ...)
    } # Fine del loop 'for (k ...)'
    
    # ... (codice invariato per il conteggio dei residui) ...
    # Questo blocco è corretto e non necessita modifiche, 
    # calcola i residui totali *presenti* in ogni ring (1, 2, 3)
    for (step in 1:Nsteps) {
      r_start <- borders[step]
      r_end <- borders[step + 1]
      selected_residues <- names(vet_dist_centroid[vet_dist_centroid >= r_start & vet_dist_centroid < r_end])
      centroids_df$residue_key <- paste(centroids_df$resid, centroids_df$resno, centroids_df$chain, sep = "_")
      df_selected <- centroids_df[centroids_df$residue_key %in% selected_residues, ]
      
      df_selected_antibody_df <- df_selected[df_selected$chain %in% chain_HL, ]
      df_selected_ligand_df <- df_selected[df_selected$chain %in% chain_AG, ]
      
      local_ab_counts[[step]] <- local_ab_counts[[step]] + table(factor(df_selected_antibody_df$resid, levels = amino_acids))
      local_ag_counts[[step]] <- local_ag_counts[[step]] + table(factor(df_selected_ligand_df$resid, levels = amino_acids))
    }
    
    # Il return ora contiene la lista nominata di 6 matrici
    return(list(
      matrices_asym = contact_matrices_asym_rings, 
      matrices_sym = contact_matrices_sym_rings,
      ab_counts = local_ab_counts, # Lista di 3 vettori
      ag_counts = local_ag_counts  # Lista di 3 vettori
    ))
    
  }, error = function(e) {
    cat("Error in", pdbL[s], ":", e$message, "\n", file = "errors.log", append = TRUE)
  })
}

# Execute (CORRETTO)
for (s in seq_along(pdbL)) {
  result <- process_pdb_file(s)
  if (is.null(result)) next # Skip file if it returned NULL
  
  # NUOVO: Loop corretto per le 6 matrici di contatto
  # Cicliamo sui nomi ("1-1", "1-2", ecc.)
  for (pair in ring_pairs) { 
    contact_matrices_asym_rings_sum[[pair]] <- contact_matrices_asym_rings_sum[[pair]] + result$matrices_asym[[pair]]
    contact_matrices_sym_rings_sum[[pair]] <- contact_matrices_sym_rings_sum[[pair]] + result$matrices_sym[[pair]]
  }
  
  # Loop invariato (e corretto) per i 3 vettori di conteggio
  for (step in 1:Nsteps) {
    residue_counts_interface_ab_sum[[step]] <- residue_counts_interface_ab_sum[[step]] + result$ab_counts[[step]]
    residue_counts_interface_ligand_df_sum[[step]] <- residue_counts_interface_ligand_df_sum[[step]] + result$ag_counts[[step]]
  }
}

# Save results (CORRETTO)
datasets <- list(
  contact_matrices_asym_rings = contact_matrices_asym_rings_sum,
  contact_matrices_sym_rings = contact_matrices_sym_rings_sum,
  residue_counts_interface_ab = residue_counts_interface_ab_sum,
  residue_counts_interface_ligand_df = residue_counts_interface_ligand_df_sum
)

# MODIFICATO: Loop di salvataggio
for (name in names(datasets)) {
  
  # Salva l'intero oggetto lista (questo era già corretto)
  saveRDS(datasets[[name]], file = file.path(results_dir_intra, paste0(name, ".rds")))
  
  # Distingui come salvare i CSV
  if (name %in% c("contact_matrices_asym_rings", "contact_matrices_sym_rings")) {
    # Caso 1: Le matrici (lista di 6 elementi nominati)
    for (pair in ring_pairs) {
      # Usiamo il nome della coppia (es. "1-1", "1-2") nel nome del file
      fname <- file.path(results_dir_intra, paste0(name, "_", pair, ".csv"))
      write.csv(datasets[[name]][[pair]], file = fname, row.names = TRUE)
    }
  } else {
    # Caso 2: I conteggi (lista di 3 elementi numerati)
    for (step in 1:Nsteps) {
      fname <- file.path(results_dir_intra, paste0(name, "_step", step, ".csv"))
      write.csv(datasets[[name]][[step]], file = fname, row.names = TRUE)
    }
  }
}

### CHECK ON THE CONTACT NUMBERS IN EACH RADIAL RING ###

### Load saved data ###
contact_matrix_sym <- readRDS(file.path(results_dir_intra, "contact_matrices_sym_rings.rds"))
contact_matrix_asym <- readRDS(file.path(results_dir_intra, "contact_matrices_asym_rings.rds"))
residue_counts_interface_ab <- readRDS(file.path(results_dir_intra, "residue_counts_interface_ab.rds"))
residue_counts_interface_ligando <- readRDS(file.path(results_dir_intra, "residue_counts_interface_ligand_df.rds"))

### Calculate and display triangular sums of contact matrices ###
triangular_sums <- numeric(2*Nsteps)

for (step in 1:(2*Nsteps)) {
  mat <- contact_matrix_sym[[step]]
  lower_tri <- mat[lower.tri(mat, diag = TRUE)]
  triangular_sums[step] <- sum(lower_tri)
  cat("Lower triangular sum for ring", step, ":", triangular_sums[step], "\n")
  
  asym_mat <- contact_matrix_asym[[step]]
  asym_mat_sum <- sum(asym_mat)
  cat("Asym mat sum for ring", step, ":", triangular_sums[step], "\n")
  
  
}

total_triangular_sum <- sum(triangular_sums) #In each ring, a number of contact >= than 44.000 (in our case) must be present to ensure statistical significativity

#Lower triangular sum for ring 1 : 44776 
#Asym mat sum for ring 1 : 44776 
#Lower triangular sum for ring 2 : 80754 
#Asym mat sum for ring 2 : 80754 
#Lower triangular sum for ring 3 : 45398 
#Asym mat sum for ring 3 : 45398 
#Lower triangular sum for ring 4 : 86994 
#Asym mat sum for ring 4 : 86994 
#Lower triangular sum for ring 5 : 70882 
#Asym mat sum for ring 5 : 70882 
#Lower triangular sum for ring 6 : 9342 
#Asym mat sum for ring 6 : 9342 