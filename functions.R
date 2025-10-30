### RENUMBERING FUNCTION ###

# Richiede: bio3d
# Input:  pdb_aus  -> oggetto ritornato da read.pdb()
#         pdb_path -> percorso del file (solo per logging/messaggi)
# Output: list(ok=TRUE, df_coord_renumbered=...)  oppure ok=FALSE + error

renumber_ab_chains <- function(pdb_aus, pdb_path, log_file = "errors.log") {
  
  # To speed up contact calculations, we use a coarse-grained representation:
  # we retain only the centroids of side chains for each residue.
  # This excludes backbone atoms (N, CA, C, O) for all amino acids,
  # except glycine, where we keep the CA atom (since its side chain is just a hydrogen).
  df_coord <- pdb_aus$atom[pdb_aus$atom$type != "HETATM" & (
    (pdb_aus$atom$resid != "GLY" & !(pdb_aus$atom$elety %in% c("N", "CA", "C", "O"))) |
      (pdb_aus$atom$resid == "GLY" & pdb_aus$atom$elety == "CA")),]
  
  if (nrow(df_coord) == 0) {
    writeLines(paste("File whose df_coord is empty:", pdb_path), log_file, append = TRUE)
    return(list(ok = FALSE, error = "df_coord empty"))
  }
  
  # Create a copy of the coordinates to reassign residue numbers per chain
  # This step resolves issues due to non-standard residue numbering in antibodies, e.g. insertions like SER100, LYS100A, THR100B, ASN101...
  df_coord_renumbered <- df_coord
  corrected_dfs <- list()
  
  # Reassign residue numbers sequentially within each chain
  for (chain in unique(df_coord_renumbered$chain)) {
    residui_unici <- unique(paste(df_coord_renumbered$resno[df_coord_renumbered$chain == chain], 
                                  df_coord_renumbered$insert[df_coord_renumbered$chain == chain], sep=""))
    nuova_numerazione <- setNames(seq_along(residui_unici), residui_unici)
    df_coord_renumbered$resno[df_coord_renumbered$chain == chain] <- nuova_numerazione[paste(df_coord_renumbered$resno[df_coord_renumbered$chain == chain], 
                                                                                           df_coord_renumbered$insert[df_coord_renumbered$chain == chain], sep="")]
    corrected_dfs[[chain]] <- df_coord_renumbered[df_coord_renumbered$chain == chain, ]
  }
  
  df_coord_renumbered <- do.call(rbind, corrected_dfs) # Merge corrected data from all chains into a single dataframe
  
  return(list(ok = TRUE, df_coord_renumbered = df_coord_renumbered))
}

### COMPUTING AVERAGE HYDROPATHY COMPLEMENTARITY ###

#Function to process each PDB file in the pdbL list and calculate an average hydropathy complementarity value for the whole interface
avg_hydropathy_complementarity <- function(pdb_path) {
  
  file_name <- basename(pdb_path)
  
  tryCatch({
    cat("Elaborating:", file_name, "\n")
    
    # Retrieve and read the PDB file
    pdb_aus <- read.pdb(pdb_path)
    
    # The following step resolves issues due to non-standard residue numbering in antibodies, e.g. insertions like SER100, LYS100A, THR100B, ASN101...
    renumbered_df <- renumber_ab_chains(pdb_aus, pdb_path, log_file = "errors.log")
    if (!renumbered_df$ok) {
      return(list(ok = FALSE, filename = file_name, path = pdb_path, error = renumbered_df$error))
    }
    df_coord_renumbered <- renumbered_df$df_coord_renumbered
    
    # Compute centroid coordinates for each residue (average of all atom coordinates)
    centroidi_df <- as.data.frame(
      df_coord_renumbered %>%
        group_by(chain, resno, resid) %>%
        summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
    )
    
    res_names <- paste(centroidi_df$resid, centroidi_df$resno, centroidi_df$chain, sep = "_")
    df_coord_resid_xyz <- centroidi_df[, c("resid", "resno", "x", "y", "z")]
    rownames(df_coord_resid_xyz) <- res_names
    
    DistMat <- as.matrix(dist(df_coord_resid_xyz[, c("x", "y", "z")])) # Compute the pairwise Euclidean distance matrix between residue centroids
    DistMat_bin <- ifelse(DistMat <= DistCutoff, 1, 0) # Binarize the distance matrix: 1 if distance ≤ cutoff, 0 otherwise
    
    # Define logical conditions to separate antigen and antibody residues
    condA <- !(centroidi_df$chain %in% c("H", "L")) # Antigen residues
    condHL <- centroidi_df$chain %in% c("H", "L") # Antibody residues
    
    # Subset the distance matrix to contain only antigen–antibody distances
    Inter_DistMat <- DistMat[condA, condHL]
    Inter_DistMat_Bin <- ifelse(Inter_DistMat <= DistCutoff, 1, 0)
    
    BS_A <- apply(Inter_DistMat_Bin, 1, sum)
    BS_A <- BS_A[BS_A != 0]
    BS_HL <- apply(Inter_DistMat_Bin, 2, sum)
    BS_HL <- BS_HL[BS_HL != 0]
    
    # Get names of residues involved in the interface
    BS_A_names <- names(BS_A)
    BS_HL_names <- names(BS_HL)
    
    # Initialize an amino acid – amino acid contact matrix (20×20)
    contact_matrix <- matrix(0, nrow=20, ncol=20, dimnames=list(amino_acids, amino_acids))
    
    # Loop through each antibody residue in contact
    # By convention, each contact is counted as bidirectional (i.e., counted twice: once for each residue involved)
    for (res_ab in BS_HL_names) {
      idx_ab <- which(colnames(Inter_DistMat_Bin) == res_ab)
      contacts <- which(Inter_DistMat_Bin[, idx_ab] == 1)
      if (length(contacts) > 0) {
        for (idx_ag in contacts) {
          res_ag <- rownames(Inter_DistMat_Bin)[idx_ag]
          aa_ab <- strsplit(res_ab, "_")[[1]][1]
          aa_ag <- strsplit(res_ag, "_")[[1]][1]
          if (aa_ab %in% amino_acids && aa_ag %in% amino_acids) {
            if(aa_ab == aa_ag){
              contact_matrix[aa_ab, aa_ag] <- contact_matrix[aa_ab, aa_ag] + 2 # Increase diagonal entry for homotypic contact
            } else{
              # Update both symmetric entries for heterotypic contact
              contact_matrix[aa_ab, aa_ag] <- contact_matrix[aa_ab, aa_ag] + 2
              contact_matrix[aa_ag, aa_ab] <- contact_matrix[aa_ag, aa_ab] + 2
            }
          }
        }
      }
    }
    
    # Create a lower-triangular version of the contact matrix for further analysis
    emi_matrix <- contact_matrix
    emi_matrix[!lower.tri(contact_matrix, diag = TRUE)] <- 0
    total_contacts <- sum(emi_matrix) # Total number of contacts in the interface
    
    # Apply hydropathy complementarity weights to the contact matrix to get an hydropathy score
    Hydropathy_complementarity_matrix <- contact_matrix * Hr_matrix
    triangular_matrix <- Hydropathy_complementarity_matrix
    triangular_matrix[!lower.tri(Hydropathy_complementarity_matrix, diag = TRUE)] <- 0
    total_hydropathy <- sum(triangular_matrix)
    avg_interface_Hr <- total_hydropathy/total_contacts
    
    return(list(
      ok = TRUE,
      filename = file_name,
      avg_interface_Hr = avg_interface_Hr
    ))
    
  }, error = function(e) {
    cat("Error in file:", pdb_path, "\n", e$message, "\n")
    return(list(ok = FALSE, filename = file_name, path = pdb_path, error = e$message))
  })
}

#Function to process each PDB file in the pdbL list and calculate an average hydropathy complementarity value for the whole interface
avg_hydropathy_complementarity_dimers <- function(pdb_path) {
  
  file_name <- basename(pdb_path)
  
  tryCatch({
    cat("Elaborating:", file_name, "\n")
    
    # Retrieve and read the PDB file
    pdb_aus <- read.pdb(pdb_path)
    
    # The following step resolves issues due to non-standard residue numbering in antibodies, e.g. insertions like SER100, LYS100A, THR100B, ASN101...
    renumbered_df <- renumber_ab_chains(pdb_aus, pdb_path, log_file = "errors.log")
    if (!renumbered_df$ok) {
      return(list(ok = FALSE, filename = file_name, path = pdb_path, error = renumbered_df$error))
    }
    df_coord_renumbered <- renumbered_df$df_coord_renumbered
    
    # Compute centroid coordinates for each residue (average of all atom coordinates)
    centroidi_df <- as.data.frame(
      df_coord_renumbered %>%
        group_by(chain, resno, resid) %>%
        summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
    )
    
    res_names <- paste(centroidi_df$resid, centroidi_df$resno, centroidi_df$chain, sep = "_")
    df_coord_resid_xyz <- centroidi_df[, c("resid", "resno", "x", "y", "z")]
    rownames(df_coord_resid_xyz) <- res_names
    
    DistMat <- as.matrix(dist(df_coord_resid_xyz[, c("x", "y", "z")])) # Compute the pairwise Euclidean distance matrix between residue centroids
    DistMat_bin <- ifelse(DistMat <= DistCutoff, 1, 0) # Binarize the distance matrix: 1 if distance ≤ cutoff, 0 otherwise
    
    # Define logical conditions to separate antigen and antibody residues
    condA <- !(centroidi_df$chain %in% "B") # Antigen residues
    condHL <- centroidi_df$chain %in% "B" # Antibody residues
    
    # Subset the distance matrix to contain only antigen–antibody distances
    Inter_DistMat <- DistMat[condA, condHL]
    Inter_DistMat_Bin <- ifelse(Inter_DistMat <= DistCutoff, 1, 0)
    
    BS_A <- apply(Inter_DistMat_Bin, 1, sum)
    BS_A <- BS_A[BS_A != 0]
    BS_HL <- apply(Inter_DistMat_Bin, 2, sum)
    BS_HL <- BS_HL[BS_HL != 0]
    
    # Get names of residues involved in the interface
    BS_A_names <- names(BS_A)
    BS_HL_names <- names(BS_HL)
    
    # Initialize an amino acid – amino acid contact matrix (20×20)
    contact_matrix <- matrix(0, nrow=20, ncol=20, dimnames=list(amino_acids, amino_acids))
    
    # Loop through each antibody residue in contact
    # By convention, each contact is counted as bidirectional (i.e., counted twice: once for each residue involved)
    for (res_ab in BS_HL_names) {
      idx_ab <- which(colnames(Inter_DistMat_Bin) == res_ab)
      contacts <- which(Inter_DistMat_Bin[, idx_ab] == 1)
      if (length(contacts) > 0) {
        for (idx_ag in contacts) {
          res_ag <- rownames(Inter_DistMat_Bin)[idx_ag]
          aa_ab <- strsplit(res_ab, "_")[[1]][1]
          aa_ag <- strsplit(res_ag, "_")[[1]][1]
          if (aa_ab %in% amino_acids && aa_ag %in% amino_acids) {
            if(aa_ab == aa_ag){
              contact_matrix[aa_ab, aa_ag] <- contact_matrix[aa_ab, aa_ag] + 2 # Increase diagonal entry for homotypic contact
            } else{
              # Update both symmetric entries for heterotypic contact
              contact_matrix[aa_ab, aa_ag] <- contact_matrix[aa_ab, aa_ag] + 2
              contact_matrix[aa_ag, aa_ab] <- contact_matrix[aa_ag, aa_ab] + 2
            }
          }
        }
      }
    }
    
    # Create a lower-triangular version of the contact matrix for further analysis
    emi_matrix <- contact_matrix
    emi_matrix[!lower.tri(contact_matrix, diag = TRUE)] <- 0
    total_contacts <- sum(emi_matrix) # Total number of contacts in the interface
    
    # Apply hydropathy complementarity weights to the contact matrix to get an hydropathy score
    Hydropathy_complementarity_matrix <- contact_matrix * Hr_matrix
    triangular_matrix <- Hydropathy_complementarity_matrix
    triangular_matrix[!lower.tri(Hydropathy_complementarity_matrix, diag = TRUE)] <- 0
    total_hydropathy <- sum(triangular_matrix)
    avg_interface_Hr <- total_hydropathy/total_contacts
    
    return(list(
      ok = TRUE,
      filename = file_name,
      avg_interface_Hr = avg_interface_Hr
    ))
    
  }, error = function(e) {
    cat("Error in file:", pdb_path, "\n", e$message, "\n")
    return(list(ok = FALSE, filename = file_name, path = pdb_path, error = e$message))
  })
}

### COMPUTE HYDROPATHY COMPLEMENTARITY VALUE BETWEEN AMINO ACIDS
### Compute the hydrophobicity complementarity 20x20 matrix between each possible pair of residues ###
#### The Milanetti scale was used to retrieve hydropathy values for each amino acid ### 
# Args:
#   milanetti_scale_path : path al file MilanettiScale.csv
#   amino_acids    : vettore dei 20 aa nell'ordine desiderato
# Return:
#   list(Hr_vector = ..., Hr_matrix = ...)
hydropathy_complementarity_matrix <- function(milanetti_scale_path, amino_acids) {
  Hr_scale <- read.csv(milanetti_scale_path)
  Hr_scale$Aa_upper <- toupper(Hr_scale$Aa)
  Hr_scale$Aa_upper <- factor(Hr_scale$Aa_upper, levels = amino_acids, ordered = TRUE)
  Hr_scale_ordered <- Hr_scale[order(Hr_scale$Aa_upper), ]
  Hr_vector <- setNames(Hr_scale_ordered$Hr, Hr_scale_ordered$Aa_upper)
  
  # Function to calculate the parabolic function hydropathy complementarity as defined in Grassmann et al. 2025 https://doi.org/10.1021/acs.jcim.4c02286
  H_formula <- function(x) {
    -0.033 * (x^2) + 0.363 * x
  }
  
  product_matrix <- outer(Hr_vector, Hr_vector, FUN = "*")
  Hr_matrix <- H_formula(product_matrix)
  rownames(Hr_matrix) <- names(Hr_vector)
  colnames(Hr_matrix) <- names(Hr_vector)
  
  list(Hr_vector = Hr_vector, Hr_matrix = Hr_matrix)
}

# Function to compute basic graph descriptors
Local_Network_Des_Bin <- function(g){
  BetweennessCentrality <- betweenness(g, v= V(g), directed = FALSE, normalized = TRUE)
  ClosenessCentrality <- closeness(g, v= V(g), normalized = TRUE)
  ShortestPath <- mean_distance(g)
  Degree <- degree(g)
  ClusteringCoefficient <- transitivity(g, type = "global")
  Density <- edge_density(g, loops = FALSE)
  Modularity <- modularity(g, membership(cluster_walktrap(g)))
  
  df_des <- as.data.frame(cbind(BetweennessCentrality, ClosenessCentrality, ShortestPath, Degree, ClusteringCoefficient, Density, Modularity))
  return((df_des))
}

# Function to compute graph energy and adjacency matrix based on spatial proximity
compute_graph_energy <- function(df, cutoff) {
  dist_mat <- as.matrix(dist(df[, c("x", "y", "z")]))
  adj <- ifelse(dist_mat <= cutoff, 1, 0)
  diag(adj) <- 0
  eig <- eigen(adj)
  list(energy = sum(abs(eig$values)), adj = adj)
}

#Function to process each PDB file in the pdbL list
geometrical_topological_descriptors_old <- function(pdb_path) {
  
  file_name <- basename(pdb_path)
  
  tryCatch({
    cat("Elaborating:", file_name, "\n")
    
    # Retrieve and read the PDB file
    pdb_aus <- read.pdb(pdb_path)
    
    # The following step resolves issues due to non-standard residue numbering in antibodies, e.g. insertions like SER100, LYS100A, THR100B, ASN101...
    renumbered_df <- renumber_ab_chains(pdb_aus, pdb_path, log_file = "errors.log")
    if (!renumbered_df$ok) {
      return(list(ok = FALSE, filename = file_name, path = pdb_path, error = renumbered_df$error))
    }
    df_coord_renumbered <- renumbered_df$df_coord_renumbered
    
    # Compute centroid coordinates for each residue (average of all atom coordinates)
    centroid_df <- as.data.frame(
      df_coord_renumbered %>%
        group_by(chain, resno, resid) %>%
        summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
    )
    
    res_names <- paste(centroid_df$resid, centroid_df$resno, centroid_df$chain, sep = "_")
    df_coord_resid_xyz <- centroid_df[, c("resid", "resno", "x", "y", "z")]
    rownames(df_coord_resid_xyz) <- res_names
    
    # Define logical conditions to separate antigen and antibody residues
    condA <- !(centroid_df$chain %in% c("H", "L")) # Antigen residues
    condHL <- centroid_df$chain %in% c("H", "L") # Antibody residues
    
    DistMat <- as.matrix(dist(df_coord_resid_xyz[, c("x", "y", "z")])) # Compute the pairwise Euclidean distance matrix between residue centroids
    DistMat_bin <- ifelse(DistMat <= DistCutoff, 1, 0) # Binarize the distance matrix: 1 if distance ≤ cutoff, 0 otherwise
    
    # Subset the distance matrix to contain only antigen–antibody distances
    Inter_DistMat <- DistMat[condA, condHL]
    Inter_DistMat_Bin <- ifelse(Inter_DistMat <= DistCutoff, 1, 0)
    
    BS_A <- apply(Inter_DistMat_Bin, 1, sum)
    BS_A <- BS_A[BS_A != 0]
    BS_HL <- apply(Inter_DistMat_Bin, 2, sum)
    BS_HL <- BS_HL[BS_HL != 0]
    
    # Get names of residues involved in the interface
    BS_A_names <- names(BS_A)
    BS_HL_names <- names(BS_HL)
    
    # --- Topological descriptors ---
    
    interface_residues <- union(BS_A_names, BS_HL_names)
    interface_centroid_df <- df_coord_resid_xyz[rownames(df_coord_resid_xyz) %in% interface_residues, ]
    
    # Compute graph energy values for interface, antigen, and antibody
    graph_energy <- compute_graph_energy(interface_centroid_df, DistCutoff)
    antigen_graph_energy <- compute_graph_energy(df_coord_resid_xyz[rownames(df_coord_resid_xyz) %in% BS_A_names, ], DistCutoff)
    ab_graph_energy <- compute_graph_energy(df_coord_resid_xyz[rownames(df_coord_resid_xyz) %in% BS_HL_names, ], DistCutoff)
    
    # Compute complex formation energy
    graph_energy_of_complexes <- graph_energy$energy - (antigen_graph_energy$energy + ab_graph_energy$energy)
    
    # Create connected, undirected graph and remove isolated nodes
    network <- graph_from_adjacency_matrix(graph_energy$adj, mode = "undirected") 
    Isolated = which(degree(network)==0)
    network_all_connected = delete_vertices(network, Isolated)
    
    # Compute network-based descriptors through the previously defined function
    df_Net <- Local_Network_Des_Bin(network_all_connected)
    df_Net$ClosenessCentrality[!is.finite(df_Net$ClosenessCentrality)] <- 0
    
    # Remove nodes with zero degree (unconnected)
    df_Net_BS <- df_Net[df_Net$Degree != 0, ]
    col_aus <- colnames(df_Net_BS)
    mean_NetDes <- apply(df_Net_BS[,col_aus], 2, mean) # Average network descriptors
    mean_NetDes_df <- as.data.frame(t(mean_NetDes))
    
    # Aggiungi le nuove metriche al dataframe
    mean_NetDes_df$GraphEnergy <- graph_energy$energy
    mean_NetDes_df$GraphEnergyOfComplexes <- graph_energy_of_complexes
    mean_NetDes_df <- cbind(PDB_file = file_name, mean_NetDes_df)
    
    # --- Geometric and numeric descriptors ---
    
    df_centroids_BS <- df_coord_resid_xyz[rownames(df_coord_resid_xyz) %in% c(BS_A_names, BS_HL_names), ]
    center_BS <- apply(df_centroids_BS[, c("x", "y", "z")], 2, mean)
    
    DistMat_centroid <- as.matrix(dist(rbind(df_centroids_BS[, c("x", "y", "z")], center_BS)))
    vet_dist_centroid <- DistMat_centroid[-nrow(DistMat_centroid), ncol(DistMat_centroid)]
    
    if (length(vet_dist_centroid) == 0) {
      writeLines(paste("No contacts found for file:", pdb_path), "errors.log", append = TRUE)
      return(list(ok = FALSE, filename = file_name, path = pdb_path, error = "No residues within cutoff distance"))
    }
    
    antigen_length <- length(unique(centroid_df$resno[centroid_df$chain == "A"])) # Differentiate between peptide-binding Ag and protein-binding Ag
    interface_radius <- round(max(vet_dist_centroid), 3)
    paratope_aa_number <- length(BS_HL)
    epitope_aa_number <- length(BS_A)
    #interface_gyration_radius <- round(sqrt(mean(rowSums((df_centroids_BS[, c("x", "y", "z")] - center_BS)^2))), 2)
    interface_gyration_radius <- round(rgyr(as.matrix(df_centroids_BS[, c("x","y","z")])), 3)
    
    
    # Add descriptors to output
    mean_NetDes_df$antigen_length <- antigen_length
    mean_NetDes_df$interface_radius <- interface_radius
    mean_NetDes_df$paratope_aa_number <- paratope_aa_number
    mean_NetDes_df$epitope_aa_number <- epitope_aa_number
    mean_NetDes_df$interface_gyration_radius <- interface_gyration_radius
    
    return(list(
      ok = TRUE,
      filename = file_name,
      mean_NetDes = mean_NetDes_df
    ))
    
  }, error = function(e) {
    cat("Error in file:", pdb_path, "\n", e$message, "\n")
    return(list(ok = FALSE, filename = file_name, path = pdb_path, error = e$message))
  })
}

#Function to process each PDB file in the pdbL list
geometrical_topological_descriptors <- function(pdb_path) {
  
  file_name <- basename(pdb_path)
  
  tryCatch({
    cat("Elaborating:", file_name, "\n")
    
    # Retrieve and read the PDB file
    pdb_aus <- read.pdb(pdb_path)
    
    # The following step resolves issues due to non-standard residue numbering in antibodies, e.g. insertions like SER100, LYS100A, THR100B, ASN101...
    renumbered_df <- renumber_ab_chains(pdb_aus, pdb_path, log_file = "errors.log")
    if (!renumbered_df$ok) {
      return(list(ok = FALSE, filename = file_name, path = pdb_path, error = renumbered_df$error))
    }
    df_coord_renumbered <- renumbered_df$df_coord_renumbered
    
    # Compute centroid coordinates for each residue (average of all atom coordinates)
    centroid_df <- as.data.frame(
      df_coord_renumbered %>%
        group_by(chain, resno, resid) %>%
        summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
    )
    
    res_names <- paste(centroid_df$resid, centroid_df$resno, centroid_df$chain, sep = "_")
    df_coord_resid_xyz <- centroid_df[, c("resid", "resno", "x", "y", "z")]
    rownames(df_coord_resid_xyz) <- res_names
    
    # Define logical conditions to separate antigen and antibody residues
    condA <- !(centroid_df$chain %in% c("H", "L")) # Antigen residues
    condHL <- centroid_df$chain %in% c("H", "L") # Antibody residues
    
    DistMat <- as.matrix(dist(df_coord_resid_xyz[, c("x", "y", "z")])) # Compute the pairwise Euclidean distance matrix between residue centroids
    DistMat_bin <- ifelse(DistMat <= DistCutoff, 1, 0) # Binarize the distance matrix: 1 if distance ≤ cutoff, 0 otherwise
    
    # Subset the distance matrix to contain only antigen–antibody distances
    Inter_DistMat <- DistMat[condA, condHL]
    Inter_DistMat_Bin <- ifelse(Inter_DistMat <= DistCutoff, 1, 0)
    
    BS_A <- apply(Inter_DistMat_Bin, 1, sum)
    BS_A <- BS_A[BS_A != 0]
    BS_HL <- apply(Inter_DistMat_Bin, 2, sum)
    BS_HL <- BS_HL[BS_HL != 0]
    
    # Get names of residues involved in the interface
    BS_A_names <- names(BS_A)
    BS_HL_names <- names(BS_HL)
    
    # --- Topological descriptors ---
    
    interface_residues <- union(BS_A_names, BS_HL_names)
    interface_centroid_df <- df_coord_resid_xyz[rownames(df_coord_resid_xyz) %in% interface_residues, ]
    
    # Compute graph energy values for interface, antigen, and antibody
    graph_energy <- compute_graph_energy(interface_centroid_df, DistCutoff)
    antigen_graph_energy <- compute_graph_energy(df_coord_resid_xyz[rownames(df_coord_resid_xyz) %in% BS_A_names, ], DistCutoff)
    ab_graph_energy <- compute_graph_energy(df_coord_resid_xyz[rownames(df_coord_resid_xyz) %in% BS_HL_names, ], DistCutoff)
    
    # Compute complex formation energy
    graph_energy_of_complexes <- graph_energy$energy - (antigen_graph_energy$energy + ab_graph_energy$energy)
    
    # Create connected, undirected graph and remove isolated nodes
    network <- graph_from_adjacency_matrix(graph_energy$adj, mode = "undirected") 
    Isolated = which(degree(network)==0)
    network_all_connected = delete_vertices(network, Isolated)
    
    # --- Calcolo diretto dei descrittori (inlining della ex Local_Network_Des_Bin) ---
    betweenness_centrality <- betweenness(network_all_connected, v = V(network_all_connected), directed = FALSE, normalized = TRUE)
    closeness_centrality <- closeness(network_all_connected, v = V(network_all_connected), normalized = TRUE)
    shortest_path <- mean_distance(network_all_connected)
    degree <- degree(network_all_connected)
    clustering_coefficient <- transitivity(network_all_connected, type = "global")
    density <- edge_density(network_all_connected, loops = FALSE)
    modularity <- modularity(network_all_connected, membership(cluster_walktrap(network_all_connected)))
    
    net_des_df <- as.data.frame(cbind(
      betweenness_centrality,
      closeness_centrality,
      shortest_path,
      degree,
      clustering_coefficient,
      density,
      modularity
    ))
    
    net_des_df$closeness_centrality[!is.finite(net_des_df$closeness_centrality)] <- 0
    
    # Remove nodes with zero degree (unconnected)
    net_des_BS <- net_des_df[net_des_df$degree != 0, ]
    col_aus <- colnames(net_des_BS)
    avg_net_des <- apply(net_des_BS[,col_aus], 2, mean) # Average network descriptors
    avg_net_des_df <- as.data.frame(t(avg_net_des))
    
    # Aggiungi le nuove metriche al dataframe
    avg_net_des_df$graph_energy <- graph_energy$energy
    avg_net_des_df$graph_energy_of_complexes <- graph_energy_of_complexes
    avg_net_des_df <- cbind(PDB_file = file_name, avg_net_des_df)
    
    # --- Geometric and numeric descriptors ---
    
    df_centroids_BS <- df_coord_resid_xyz[rownames(df_coord_resid_xyz) %in% c(BS_A_names, BS_HL_names), ]
    center_BS <- apply(df_centroids_BS[, c("x", "y", "z")], 2, mean)
    
    DistMat_centroid <- as.matrix(dist(rbind(df_centroids_BS[, c("x", "y", "z")], center_BS)))
    vet_dist_centroid <- DistMat_centroid[-nrow(DistMat_centroid), ncol(DistMat_centroid)]
    
    if (length(vet_dist_centroid) == 0) {
      writeLines(paste("No contacts found for file:", pdb_path), "errors.log", append = TRUE)
      return(list(ok = FALSE, filename = file_name, path = pdb_path, error = "No residues within cutoff distance"))
    }
    
    antigen_length <- length(unique(centroid_df$resno[centroid_df$chain == "A"])) # Differentiate between peptide-binding Ag and protein-binding Ag
    interface_radius <- round(max(vet_dist_centroid), 3)
    paratope_aa_number <- length(BS_HL)
    epitope_aa_number <- length(BS_A)
    interface_gyration_radius <- round(sqrt(mean(rowSums((df_centroids_BS[, c("x", "y", "z")] - center_BS)^2))), 2)
    #interface_gyration_radius <- round(rgyr(as.matrix(df_centroids_BS[, c("x","y","z")])), 3)
    
    
    # Add descriptors to output
    avg_net_des_df$antigen_length <- antigen_length
    avg_net_des_df$interface_radius <- interface_radius
    avg_net_des_df$paratope_aa_number <- paratope_aa_number
    avg_net_des_df$epitope_aa_number <- epitope_aa_number
    avg_net_des_df$interface_gyration_radius <- interface_gyration_radius
    
    return(list(
      ok = TRUE,
      filename = file_name,
      avg_net_des_df = avg_net_des_df
    ))
    
  }, error = function(e) {
    cat("Error in file:", pdb_path, "\n", e$message, "\n")
    return(list(ok = FALSE, filename = file_name, path = pdb_path, error = e$message))
  })
}

################################################################################
################################################################################
################################################################################
# CAPRI DESCRIPTORS
################################################################################
################################################################################
################################################################################

################################################################################
# I_RMSD
################################################################################

# Restituisce un data frame con le coordinate dei residui d’interfaccia
# (righe = residui di A o H/L che hanno almeno un contatto <= DistCutoff)
find_interface_coordinates <- function(pdb_path, DistCutoff, chains_ab = c("H","L"), log_file = "errors.log") {
  file_name <- basename(pdb_path)
  tryCatch({
    cat("Elaborating:", file_name, "\n")
    
    # 1) Leggi PDB
    pdb_aus <- read.pdb(pdb_path)
    
    # 2) Rinumerazione catene Ab (H/L) se necessario
    renumbered_df <- renumber_ab_chains(pdb_aus, pdb_path, log_file = log_file)
    if (!isTRUE(renumbered_df$ok)) {
      warning(sprintf("Renumbering failed for %s: %s", file_name, renumbered_df$error))
      return(data.frame(resid=character(), resno=integer(), x=double(), y=double(), z=double()))
    }
    df_coord_renumbered <- renumbered_df$df_coord_renumbered
    
    # 3) Centroidi per residuo
    centroid_df <- as.data.frame(
      df_coord_renumbered %>%
        dplyr::group_by(chain, resno, resid) %>%
        dplyr::summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
    )
    
    # 4) Data frame con solo (resid, resno, x, y, z) e rownames "RES_RESNO_CHAIN"
    res_names <- paste(centroid_df$resid, centroid_df$resno, centroid_df$chain, sep = "_")
    df_coord_resid_xyz <- centroid_df[, c("resid", "resno", "x", "y", "z")]
    rownames(df_coord_resid_xyz) <- res_names
    
    # 5) Maschere: antigene (non H/L) vs anticorpo (H/L)
    condA  <- !(centroid_df$chain %in% chains_ab)  # Antigen
    condHL <-  (centroid_df$chain %in% chains_ab)  # Antibody
    
    if (!any(condA) || !any(condHL)) {
      warning(sprintf("No antigen or antibody residues found in %s.", file_name))
      return(df_coord_resid_xyz[0, , drop = FALSE])
    }
    
    # 6) Matrice delle distanze tra centroidi
    DistMat <- as.matrix(dist(df_coord_resid_xyz[, c("x", "y", "z")]))
    
    # 7) Solo distanze Ag–Ab e binarizzazione con cutoff
    Inter_DistMat      <- DistMat[condA, condHL, drop = FALSE]
    if (length(Inter_DistMat) == 0) {
      return(df_coord_resid_xyz[0, , drop = FALSE])
    }
    Inter_DistMat_Bin  <- ifelse(Inter_DistMat <= DistCutoff, 1, 0)
    
    # 8) Residui dell’interfaccia = quelli con almeno un contatto
    BS_A  <- rowSums(Inter_DistMat_Bin); BS_A  <- BS_A[BS_A != 0]
    BS_HL <- colSums(Inter_DistMat_Bin); BS_HL <- BS_HL[BS_HL != 0]
    
    interface_residues <- union(names(BS_A), names(BS_HL))
    
    # 9) Subset delle sole coordinate d’interfaccia
    interface_centroid_df <- df_coord_resid_xyz[rownames(df_coord_resid_xyz) %in% interface_residues,
                                                , drop = FALSE]
    
    return(interface_centroid_df)
    
  }, error = function(e) {
    warning(sprintf("Error in file %s: %s", file_name, e$message))
    # logga su file (append)
    write(sprintf("Error in %s: %s", pdb_path, e$message), file = log_file, append = TRUE)
    # Ritorna df vuoto con le colonne attese
    return(data.frame(resid=character(), resno=integer(), x=double(), y=double(), z=double()))
  })
}

################################################################################
# F_NAT
################################################################################

#Function to process each PDB file in the pdbL list and calculate an average hydropathy complementarity value for the whole interface
compute_contacts_from_file <- function(pdb_path, DistCutoff) {
  file_name <- basename(pdb_path)
  tryCatch({
    cat("Elaborating:", file_name, "\n")
    
    # Retrieve and read the PDB file
    pdb_aus <- read.pdb(pdb_path)
    
    # To speed up contact calculations, we use a coarse-grained representation:
    # we retain only the centroids of side chains for each residue.
    # This excludes backbone atoms (N, CA, C, O) for all amino acids,
    # except glycine, where we keep the CA atom (since its side chain is just a hydrogen).
    df_coord <- pdb_aus$atom[pdb_aus$atom$type != "HETATM" & (
      (pdb_aus$atom$resid != "GLY" & !(pdb_aus$atom$elety %in% c("N", "CA", "C", "O"))) |
        (pdb_aus$atom$resid == "GLY" & pdb_aus$atom$elety == "CA")),]
    
    
    if (nrow(df_coord) == 0) {
      writeLines(paste("File whose df_coord is empty:", pdb_path), "errors.log", append = TRUE)
      return(list(ok = FALSE, filename = file_name, path = pdb_path, error = "df_coord empty"))
    }
    
    # Create a copy of the coordinates to reassign residue numbers per chain
    # This step resolves issues due to non-standard residue numbering in antibodies, e.g. insertions like SER100, LYS100A, THR100B, ASN101...
    df_coord_corrected <- df_coord
    corrected_dfs <- list()
    
    # Reassign residue numbers sequentially within each chain
    for (chain in unique(df_coord_corrected$chain)) {
      residui_unici <- unique(paste(df_coord_corrected$resno[df_coord_corrected$chain == chain], 
                                    df_coord_corrected$insert[df_coord_corrected$chain == chain], sep=""))
      nuova_numerazione <- setNames(seq_along(residui_unici), residui_unici)
      df_coord_corrected$resno[df_coord_corrected$chain == chain] <- nuova_numerazione[paste(df_coord_corrected$resno[df_coord_corrected$chain == chain], 
                                                                                             df_coord_corrected$insert[df_coord_corrected$chain == chain], sep="")]
      corrected_dfs[[chain]] <- df_coord_corrected[df_coord_corrected$chain == chain, ]
    }
    
    df_coord_corrected <- do.call(rbind, corrected_dfs) # Merge corrected data from all chains into a single dataframe
    
    # Compute centroid coordinates for each residue (average of all atom coordinates)
    centroidi_df <- as.data.frame(
      df_coord_corrected %>%
        group_by(chain, resno, resid) %>%
        summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
    )
    
    res_names <- paste(centroidi_df$resid, centroidi_df$resno, centroidi_df$chain, sep = "_")
    df_coord_resid_xyz <- centroidi_df[, c("resid", "resno", "x", "y", "z")]
    rownames(df_coord_resid_xyz) <- res_names
    
    DistMat <- as.matrix(dist(df_coord_resid_xyz[, c("x", "y", "z")])) # Compute the pairwise Euclidean distance matrix between residue centroids
    #DistMat_bin <- ifelse(DistMat <= DistCutoff, 1, 0) # Binarize the distance matrix: 1 if distance ≤ cutoff, 0 otherwise
    
    # Define logical conditions to separate antigen and antibody residues
    condA <- !(centroidi_df$chain %in% c("H", "L")) # Antigen residues
    condHL <- centroidi_df$chain %in% c("H", "L") # Antibody residues
    
    # Subset the distance matrix to contain only antigen–antibody distances
    Inter_DistMat <- DistMat[condA, condHL]
    Inter_DistMat_Bin <- ifelse(Inter_DistMat <= DistCutoff, 1, 0)
    
    BS_A <- apply(Inter_DistMat_Bin, 1, sum)
    BS_A <- BS_A[BS_A != 0]
    BS_HL <- apply(Inter_DistMat_Bin, 2, sum)
    BS_HL <- BS_HL[BS_HL != 0]
    
    # Get names of residues involved in the interface
    BS_A_names <- names(BS_A)
    BS_HL_names <- names(BS_HL)
    
    contact_indices <- which(Inter_DistMat_Bin == 1, arr.ind = TRUE)
    contacts <- data.frame(
      res1 = rownames(Inter_DistMat_Bin)[contact_indices[, 1]],
      res2 = colnames(Inter_DistMat_Bin)[contact_indices[, 2]]
    )
    
    return(list(ok = TRUE, filename = file_name, path = pdb_path, contacts = contacts))
    
  }, error = function(e) {
    cat("Error in file:", pdb_path, "\n", e$message, "\n")
    return(list(ok = FALSE, filename = file_name, path = pdb_path, error = e$message))
  })
}

compute_f_nat <- function(reference_csv, docking_csv, unordered = FALSE){
  
  ref  <- fread(reference_csv)   # already in data.table format
  dock <- fread(docking_csv)
  
  inter <- fintersect(ref, dock)
  
  n_ref    <- nrow(ref)
  n_native <- nrow(inter)
  fnc <- n_native / n_ref
  
  list(
    n_reference = n_ref,
    n_docking   = nrow(dock),
    n_native    = n_native,
    fraction    = fnc,
    fraction_percent = round(100 * fnc, 2),
    native_contacts = inter
  )
}

################################################################################
# L_RMSD
################################################################################















#Function to process each PDB file in the pdbL list and calculate an average hydropathy complementarity value for the whole interface
cont_mat_from_dist_mat <- function(pdb_path) {
  
  file_name <- basename(pdb_path)
  
  tryCatch({
    cat("Elaborating:", file_name, "\n")
    
    # Retrieve and read the PDB file
    pdb_aus <- read.pdb(pdb_path)
    
    # The following step resolves issues due to non-standard residue numbering in antibodies, e.g. insertions like SER100, LYS100A, THR100B, ASN101...
    renumbered_df <- renumber_ab_chains(pdb_aus, pdb_path, log_file = "errors.log")
    if (!renumbered_df$ok) {
      return(list(ok = FALSE, filename = file_name, path = pdb_path, error = renumbered_df$error))
    }
    df_coord_renumbered <- renumbered_df$df_coord_renumbered
    
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
    
    # Assign contacts to rings (only intra-ring contacts considered)
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
    contact_matrices_asym_rings <- vector("list", Nsteps)
    contact_matrices_sym_rings  <- vector("list", Nsteps)
    local_ab_counts <- vector("list", Nsteps)
    local_ag_counts <- vector("list", Nsteps)
    
    for (step in 1:Nsteps) {
      contact_matrices_asym_rings[[step]] <- matrix(0, nrow=length(amino_acids), ncol=length(amino_acids), dimnames=list(amino_acids, amino_acids))
      contact_matrices_sym_rings[[step]]  <- matrix(0, nrow=length(amino_acids), ncol=length(amino_acids), dimnames=list(amino_acids, amino_acids))
      local_ab_counts[[step]] <- setNames(rep(0, length(amino_acids)), amino_acids)
      local_ag_counts[[step]] <- setNames(rep(0, length(amino_acids)), amino_acids)
    }
    
    
    # Fill intra-ring contact matrices
    for (k in seq_len(nrow(contacts))) {
      aa1 <- strsplit(contacts$res1[k], "_")[[1]][1]
      aa2 <- strsplit(contacts$res2[k], "_")[[1]][1]
      if (!(aa1 %in% amino_acids && aa2 %in% amino_acids)) next
      
      ring1 <- contacts$ring_res1[k]
      ring2 <- contacts$ring_res2[k]
      
      if (ring1 == ring2 && ring1 >= 1 && ring1 <= Nsteps) {
        assigned_ring <- ring1
        
        chain1 <- strsplit(contacts$res1[k], "_")[[1]][3]
        chain2 <- strsplit(contacts$res2[k], "_")[[1]][3]
        
        ## ---- ASYMMETRIC CONTACT MATRICES ----
        
        if (chain1 %in% chain_HL && chain2 %in% chain_AG) {
          # antibody → antigen
          contact_matrices_asym_rings[[assigned_ring]][aa1, aa2] <- contact_matrices_asym_rings[[assigned_ring]][aa1, aa2] + 2
        } else if (chain1 %in% chain_AG && chain2 %in% chain_HL) {
          # antigen → antibody (switch to keep rows = Ab, columns = Ag)
          contact_matrices_asym_rings[[assigned_ring]][aa2, aa1] <- contact_matrices_asym_rings[[assigned_ring]][aa2, aa1] + 2
        }
        
        ## ---- SYMMETRIC CONTACT MATRICES ----
        
        if (aa1 == aa2) {
          contact_matrices_sym_rings[[assigned_ring]][aa1, aa1] <- contact_matrices_sym_rings[[assigned_ring]][aa1, aa1] + 2
        } else {
          contact_matrices_sym_rings[[assigned_ring]][aa1, aa2] <- contact_matrices_sym_rings[[assigned_ring]][aa1, aa2] + 2
          contact_matrices_sym_rings[[assigned_ring]][aa2, aa1] <- contact_matrices_sym_rings[[assigned_ring]][aa2, aa1] + 2
        }
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
    
    return(list(
      matrices_asym = contact_matrices_asym_rings,
      matrices_sym  = contact_matrices_sym_rings,
      ab_counts = local_ab_counts,
      ag_counts = local_ag_counts
    ))
    
  }, error = function(e) {
    cat("Error in file:", pdb_path, "\n", e$message, "\n")
    return(list(ok = FALSE, filename = file_name, path = pdb_path, error = e$message))
  })
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
    
    # Assign contacts to rings (only intra-ring contacts considered)
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
    contact_matrices_asym_rings <- vector("list", Nsteps)
    contact_matrices_sym_rings  <- vector("list", Nsteps)
    local_ab_counts <- vector("list", Nsteps)
    local_ag_counts <- vector("list", Nsteps)
    
    for (step in 1:Nsteps) {
      contact_matrices_asym_rings[[step]] <- matrix(0, nrow=length(amino_acids), ncol=length(amino_acids), dimnames=list(amino_acids, amino_acids))
      contact_matrices_sym_rings[[step]]  <- matrix(0, nrow=length(amino_acids), ncol=length(amino_acids), dimnames=list(amino_acids, amino_acids))
      local_ab_counts[[step]] <- setNames(rep(0, length(amino_acids)), amino_acids)
      local_ag_counts[[step]] <- setNames(rep(0, length(amino_acids)), amino_acids)
    }
    
    
    # Fill intra-ring contact matrices
    for (k in seq_len(nrow(contacts))) {
      aa1 <- strsplit(contacts$res1[k], "_")[[1]][1]
      aa2 <- strsplit(contacts$res2[k], "_")[[1]][1]
      if (!(aa1 %in% amino_acids && aa2 %in% amino_acids)) next
      
      ring1 <- contacts$ring_res1[k]
      ring2 <- contacts$ring_res2[k]
      
      if (ring1 == ring2 && ring1 >= 1 && ring1 <= Nsteps) {
        assigned_ring <- ring1
        
        chain1 <- strsplit(contacts$res1[k], "_")[[1]][3]
        chain2 <- strsplit(contacts$res2[k], "_")[[1]][3]
        
        ## ---- ASYMMETRIC CONTACT MATRICES ----
        
        if (chain1 %in% chain_HL && chain2 %in% chain_AG) {
          # antibody → antigen
          contact_matrices_asym_rings[[assigned_ring]][aa1, aa2] <- contact_matrices_asym_rings[[assigned_ring]][aa1, aa2] + 2
        } else if (chain1 %in% chain_AG && chain2 %in% chain_HL) {
          # antigen → antibody (switch to keep rows = Ab, columns = Ag)
          contact_matrices_asym_rings[[assigned_ring]][aa2, aa1] <- contact_matrices_asym_rings[[assigned_ring]][aa2, aa1] + 2
        }
        
        ## ---- SYMMETRIC CONTACT MATRICES ----
        
        if (aa1 == aa2) {
          contact_matrices_sym_rings[[assigned_ring]][aa1, aa1] <- contact_matrices_sym_rings[[assigned_ring]][aa1, aa1] + 2
        } else {
          contact_matrices_sym_rings[[assigned_ring]][aa1, aa2] <- contact_matrices_sym_rings[[assigned_ring]][aa1, aa2] + 2
          contact_matrices_sym_rings[[assigned_ring]][aa2, aa1] <- contact_matrices_sym_rings[[assigned_ring]][aa2, aa1] + 2
        }
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
    
    return(list(
      matrices_asym = contact_matrices_asym_rings,
      matrices_sym  = contact_matrices_sym_rings,
      ab_counts = local_ab_counts,
      ag_counts = local_ag_counts
    ))
    
    
    
  }, error = function(e) {
    cat("Error in", pdbL[s], ":", e$message, "\n", file = "errors.log", append = TRUE)
  })
}


