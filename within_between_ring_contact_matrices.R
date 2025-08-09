### THIS SCRIPT ANALYZES ANTIBODY–ANTIGEN PDB COMPLEXES, BUILDS RESIDUE–RESIDUE CONTACT MATRICES,
### COUNTS AMINO ACIDS PER CONCENTRIC RING AROUND THE INTERFACE, AND "ENRICHES" INNER RINGS:
### IF A CONTACT SPANS TWO RINGS, THE CONTACT IS CREDITED TO THE INNER RING'S COUNTS.
### THE OUTER RESIDUE OF A CONTACT SPANNING TWO RINGS IS THEN ASSIGNED TO BOTH THE INNER AND THE OUTER RING
### RESULTS ARE SAVED TO CSV/RDS PER RING, WITH A FINAL CONTACT-CONSISTENCY CHECK.

### Libraries ###
library(bio3d)
library(dplyr)

### Setting paths, directories and global variables ###
pdb_dir <- "/path/to/your/non/redundant/pdb/files/directory/"
results_dir <- "/path/to/the/directory/where/you/want/to/save/the/matrices/and/residue/count/" 
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

DistCutoff <- 8.5  # Two residues are defined as "in contact" if the centroids of their side chain are at less than 8.5 Angstroms
Nsteps <- 3 # Number of concentric "rings" around the interface

pdbL <- list.files(pdb_dir)

# Amino acids ordered according to the Kyte-Doolittle hydrophobicity scale
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE")

### MAIN: GENERATING CONTACT MATRICES AND COMPUTING RESIDUE COUNTS ###

# Initialize cumulative data structures
contact_matrices_rings_sum <- vector("list", Nsteps)
residue_counts_interface_ab_sum <- vector("list", Nsteps) # residue counts for antibody chains (H/L)
residue_counts_interface_ligand_df_sum <- vector("list", Nsteps) # residue counts for antigen chain (A)

for (step in 1:Nsteps) {
  contact_matrices_rings_sum[[step]] <- matrix(0, nrow = 20, ncol = 20, dimnames = list(amino_acids, amino_acids))
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
    # we retain only the centroids of side chains for each residue.
    # This excludes backbone atoms (N, CA, C, O) for all amino acids,
    # except glycine, where we keep the CA atom (since its side chain is just a hydrogen).
    df_coord <- pdb_aus$atom[pdb_aus$atom$type != "HETATM" & (
      (pdb_aus$atom$resid != "GLY" & !(pdb_aus$atom$elety %in% c("N", "CA", "C", "O"))) |
        (pdb_aus$atom$resid == "GLY" & pdb_aus$atom$elety == "CA")),]
    
    if (nrow(df_coord) == 0) {
      cat("File con df_coord vuoto:", path_aus, "\n", file = "errors.log", append = TRUE)
      return(list(ok = FALSE, filename = pdbL[s], error = "df_coord vuoto"))
    } 
    
    # Create a copy of the coordinates to reassign residue numbers per chain
    # This step resolves issues due to non-standard residue numbering in antibodies, e.g. insertions like SER100, LYS100A, THR100B, ASN101...
    df_coord_corrected <- df_coord
    corrected_dfs <- list()
    
    # Reassign residue numbers sequentially within each chain
    for (chain in unique(df_coord_corrected$chain)) {
      unique_residues <- unique(paste(df_coord_corrected$resno[df_coord_corrected$chain == chain], 
                                      df_coord_corrected$insert[df_coord_corrected$chain == chain], sep=""))
      new_numbering <- setNames(seq_along(unique_residues), unique_residues)
      df_coord_corrected$resno[df_coord_corrected$chain == chain] <- new_numbering[paste(df_coord_corrected$resno[df_coord_corrected$chain == chain], 
                                                                                         df_coord_corrected$insert[df_coord_corrected$chain == chain], sep="")]
      corrected_dfs[[chain]] <- df_coord_corrected[df_coord_corrected$chain == chain, ]
    }
    
    df_coord_corrected <- do.call(rbind, corrected_dfs)
    
    # Compute centroid coordinates for each residue (average of all atom coordinates)
    centroids_df <- as.data.frame(
      df_coord_corrected %>%
        group_by(chain, resno, resid) %>%
        summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
    )
    
    res_names <- paste(centroids_df$resid, centroids_df$resno, centroids_df$chain, sep = "_")
    df_coord_resid_xyz <- centroids_df[, c("resid", "resno", "x", "y", "z")]
    rownames(df_coord_resid_xyz) <- res_names
    
    DistMat <- as.matrix(dist(df_coord_resid_xyz[, c("x", "y", "z")])) # Compute the pairwise Euclidean distance matrix between residue centroids
    
    # Define logical conditions to separate antigen and antibody residues
    condA <- !(centroids_df$chain %in% chain_HL)
    condHL <- centroids_df$chain %in% chain_HL
    
    # Subset the distance matrix to contain only antigen–antibody distances
    Inter_DistMat <- DistMat[condA, condHL]
    Inter_DistMat_Bin <- ifelse(Inter_DistMat <= DistCutoff, 1, 0)
    
    BS_A <- apply(Inter_DistMat_Bin, 1, sum)
    BS_A <- BS_A[BS_A != 0]
    BS_HL <- apply(Inter_DistMat_Bin, 2, sum)
    BS_HL <- BS_HL[BS_HL != 0]
    
    if (length(BS_A) == 0 || length(BS_HL) == 0) return(NULL)
    
    # Get names of residues involved in the interface
    BS_A_names <- names(BS_A)
    BS_HL_names <- names(BS_HL)
    
    # Calculate rings around interface centroid
    df_centroids_BS <- df_coord_resid_xyz[rownames(df_coord_resid_xyz) %in% c(BS_A_names, BS_HL_names), ]
    center_BS <- colMeans(df_centroids_BS[, c("x", "y", "z")])
    DistMat_centroid <- as.matrix(dist(rbind(df_centroids_BS[, c("x", "y", "z")], center_BS)))
    vet_dist_centroid <- DistMat_centroid[-nrow(DistMat_centroid), ncol(DistMat_centroid)]
    r_max <- max(vet_dist_centroid)
    borders <- c(0, 0.36, 0.64, 1) * r_max # Ring boundaries as fractions of max radius
    
    # Assign contacts to rings
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
    
    
    # Initialize local data structures
    contact_matrices_rings <- vector("list", Nsteps)
    local_ab_counts <- vector("list", Nsteps)
    local_ag_counts <- vector("list", Nsteps)
    
    for (step in 1:Nsteps) {
      contact_matrices_rings[[step]] <- matrix(0, nrow=length(amino_acids), ncol=length(amino_acids), dimnames=list(amino_acids, amino_acids))
      local_ab_counts[[step]] <- setNames(rep(0, length(amino_acids)), amino_acids)
      local_ag_counts[[step]] <- setNames(rep(0, length(amino_acids)), amino_acids)
    }
    
    # CONTACT MATRIX FILLING (ring assignment = min(ring of the two residues))
    # Rationale: a contact bridging rings r1 and r2 is credited to the inner ring min(r1, r2).
    
    for (k in seq_len(nrow(contacts))) {
      aa1 <- strsplit(contacts$res1[k], "_")[[1]][1]
      aa2 <- strsplit(contacts$res2[k], "_")[[1]][1]
      
      if (!(aa1 %in% amino_acids && aa2 %in% amino_acids)) next
      
      ring1 <- contacts$ring_res1[k]
      ring2 <- contacts$ring_res2[k]
      assigned_ring <- min(ring1, ring2)
      if (assigned_ring < 1 || assigned_ring > Nsteps) next
      
      if (aa1 == aa2) {
        contact_matrices_rings[[assigned_ring]][aa1, aa1] <- contact_matrices_rings[[assigned_ring]][aa1, aa1] + 2
      } else {
        contact_matrices_rings[[assigned_ring]][aa1, aa2] <- contact_matrices_rings[[assigned_ring]][aa1, aa2] + 2
        contact_matrices_rings[[assigned_ring]][aa2, aa1] <- contact_matrices_rings[[assigned_ring]][aa2, aa1] + 2
      }
    }
    
    # Count residues per ring for antibody and antigen
    for (step in 1:Nsteps) {
      
      # Determine where the ring starts and where it ends
      r_start <- borders[step]
      r_end <- borders[step + 1]
      
      # Select residues whose distance from the interface centroid falls within the current ring boundaries
      selected_residues <- names(vet_dist_centroid[vet_dist_centroid >= r_start & vet_dist_centroid < r_end])
      centroids_df$residue_key <- paste(centroids_df$resid, centroids_df$resno, centroids_df$chain, sep = "_")
      df_selected <- centroids_df[centroids_df$residue_key %in% selected_residues, ]
      
      # Update amino acid counts for antibody and antigen residues in this ring
      df_selected_antibody_df <- df_selected[df_selected$chain %in% chain_HL, ]
      df_selected_ligand_df <- df_selected[df_selected$chain %in% chain_AG, ]
      
      local_ab_counts[[step]] <- local_ab_counts[[step]] + table(factor(df_selected_antibody_df$resid, levels = amino_acids))
      local_ag_counts[[step]] <- local_ag_counts[[step]] + table(factor(df_selected_ligand_df$resid, levels = amino_acids))
    }
    
    ### ENRICHING INNER RINGS WITH OUTER RESIDUES FROM CROSS-RING CONTACTS ###
    # Goal: if a contact spans two different rings (r_outer > r_inner),
    #       also CREDIT the outer residue to the COUNT of the inner ring.
    # Important: this enrichment updates ONLY the residue counts (local_ab_counts / local_ag_counts),
    #            NOT the contact matrices (which remain assigned to min(r1, r2) above).
    # Duplicate control: each residue is credited at most once per target ring.
    
    # Track which residue keys have already been credited to a given ring (to avoid duplicates)
    
    seen_residues_ab <- vector("list", Nsteps)
    seen_residues_ag <- vector("list", Nsteps)
    for (step in 1:Nsteps) {
      seen_residues_ab[[step]] <- character()
      seen_residues_ag[[step]] <- character()
    }
    
    for (k in seq_len(nrow(contacts))) {
      
      res1 <- contacts$res1[k]
      res2 <- contacts$res2[k]
      aa1 <- strsplit(res1, "_")[[1]][1]
      aa2 <- strsplit(res2, "_")[[1]][1]
      ring1 <- contacts$ring_res1[k]
      ring2 <- contacts$ring_res2[k]
      chain1 <- strsplit(res1, "_")[[1]][3]
      chain2 <- strsplit(res2, "_")[[1]][3]
      
      # Identify the outer residue (higher ring) and the inner ring to be enriched
      if (ring1 > ring2) {
        key <- res1
        aa <- aa1
        chain <- chain1
        target_ring <- ring2
      } else if (ring2 > ring1) {
        key <- res2
        aa <- aa2
        chain <- chain2
        target_ring <- ring1
      } else {
        next  # Same ring: nothing to enrich
      }
      
      # Credit the outer residue to the inner ring ONCE (per chain type and ring)
      if (aa %in% amino_acids) {
        if (chain %in% chain_HL && !(key %in% seen_residues_ab[[target_ring]])) {
          local_ab_counts[[target_ring]][aa] <- local_ab_counts[[target_ring]][aa] + 1
          seen_residues_ab[[target_ring]] <- c(seen_residues_ab[[target_ring]], key)
        } else if (chain %in% chain_AG && !(key %in% seen_residues_ag[[target_ring]])) {
          local_ag_counts[[target_ring]][aa] <- local_ag_counts[[target_ring]][aa] + 1
          seen_residues_ag[[target_ring]] <- c(seen_residues_ag[[target_ring]], key)
        }
      }
    }
    
    return(list(
      matrices = contact_matrices_rings,
      ab_counts = local_ab_counts,
      ag_counts = local_ag_counts
    ))
    
    
  }, error = function(e) {
    cat("Error in", pdbL[s], ":", e$message, "\n", file = "errors.log", append = TRUE)
  })
}

# Execute
for (s in seq_along(pdbL)) {
  result <- process_pdb_file(s)
  if (is.null(result)) next  # Skip file if it returned NULL
  
  for (step in 1:Nsteps) {
    contact_matrices_rings_sum[[step]] <- contact_matrices_rings_sum[[step]] + result$matrices[[step]]
    residue_counts_interface_ab_sum[[step]] <- residue_counts_interface_ab_sum[[step]] + result$ab_counts[[step]]
    residue_counts_interface_ligand_df_sum[[step]] <- residue_counts_interface_ligand_df_sum[[step]] + result$ag_counts[[step]]
  }
}

# Save results 
datasets <- list(
  contact_matrix_rings = contact_matrices_rings_sum,
  residue_counts_interface_ab = residue_counts_interface_ab_sum,
  residue_counts_interface_ligand_df = residue_counts_interface_ligand_df_sum
)

for (name in names(datasets)) {
  for (step in 1:Nsteps) {
    fname <- file.path(results_dir, paste0(name, "_step", step, ".csv"))
    write.csv(datasets[[name]][[step]], file = fname, row.names = TRUE)
  }
  saveRDS(datasets[[name]], file = file.path(results_dir, paste0(name, ".rds")))
}

### CHECK ON THE CONTACT NUMBERS IN EACH RADIAL RING ###

### Load saved data ###
contact_matrix <- readRDS(file.path(results_dir, "contact_matrix_rings.rds"))
residue_counts_interface_ab <- readRDS(file.path(results_dir, "residue_counts_interface_ab.rds"))
residue_counts_interface_ligando <- readRDS(file.path(results_dir, "residue_counts_interface_ligand_df.rds"))

### Calculate and display triangular sums of contact matrices ###
triangular_sums <- numeric(Nsteps)

for (step in 1:Nsteps) {
  mat <- contact_matrix[[step]]
  lower_tri <- mat[lower.tri(mat, diag = TRUE)]
  triangular_sums[step] <- sum(lower_tri)
  cat("Lower triangular sum for ring", step, ":", triangular_sums[step], "\n")
}

total_triangular_sum <- sum(triangular_sums) #In each ring, a number of contact >= than 44.000 (in our case) must be present to ensure statistical significativity
