### THIS SCRIPT EXTRACTS TOPOLOGICAL, GEOMETRIC, AND NUMERIC DESCRIPTORS FROM ANTIBODY–ANTIGEN INTERFACES BASED ON PDB STRUCTURES ###
### OUTPUT IS A SINGLE CSV FILE CONTAINING AVERAGE DESCRIPTORS PER COMPLEX ###

### Libraries ###
library(bio3d)
library(dplyr)
library(future)
library(furrr)
library(purrr)
library(progressr)
library(igraph)

### Set parallelization to speed up computations (use all available CPU but one) ###
plan(multisession, workers = parallel::detectCores() - 1)
handlers(global = TRUE)
handlers("rstudio")

### Setting paths, directories and global variables ###
pdb_dir <- "/path/to/your/non/redundant/pdb/files/directory/"
files_dir <- "/path/to/the/directory/where/you/want/to/save/the/scores/"
DistCutoff <- 8.5  # Two residues are defined as "in contact" if the centroids of their side chain are at less than 8.5 Angstroms
pdbL <- list.files(path = pdb_dir, pattern = "*.pdb")
set.seed(1234)

# Amino acids ordered according to the Kyte-Doolittle hydrophobicity scale
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE")

### Setting functions ###

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
  list(energy = sum(eig$values), adj = adj)
}

# Function to process each PDB file in the pdbL list
process_pdb_file <- function(s) {
  tryCatch({
    cat("Elaborating:", pdbL[s], "\n")
    
    # Retrieve and read the PDB file
    path_aus <- paste0(pdb_dir, pdbL[s])
    pdb_aus <- read.pdb(path_aus)
    
    # To speed up contact calculations, we use a coarse-grained representation:
    # we retain only the centroids of side chains for each residue.
    # This excludes backbone atoms (N, CA, C, O) for all amino acids,
    # except glycine, where we keep the CA atom (since its side chain is just a hydrogen).
    df_coord <- pdb_aus$atom[pdb_aus$atom$type != "HETATM" & (
      (pdb_aus$atom$resid != "GLY" & !(pdb_aus$atom$elety %in% c("N", "CA", "C", "O"))) |
        (pdb_aus$atom$resid == "GLY" & pdb_aus$atom$elety == "CA")),]
    
    if (nrow(df_coord) == 0) {
      writeLines(paste("File con df_coord vuoto:", path_aus), "errors.log", append = TRUE)
      return(list(ok = FALSE, filename = basename(path_aus), error = "df_coord vuoto"))
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
    
    df_coord_corrected <- do.call(rbind, corrected_dfs)
    
    # Compute centroid coordinates for each residue (average of all atom coordinates)
    centroid_df <- as.data.frame(
      df_coord_corrected %>%
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
    mean_NetDes_df <- cbind(PDB_file = pdbL[s], mean_NetDes_df)

    # --- Geometric and numeric descriptors ---
    
    df_centroids_BS <- df_coord_resid_xyz[rownames(df_coord_resid_xyz) %in% c(BS_A_names, BS_HL_names), ]
    center_BS <- apply(df_centroids_BS[, c("x", "y", "z")], 2, mean)
    
    DistMat_centroid <- as.matrix(dist(rbind(df_centroids_BS[, c("x", "y", "z")], center_BS)))
    vet_dist_centroid <- DistMat_centroid[-nrow(DistMat_centroid), ncol(DistMat_centroid)]
    
    if (length(vet_dist_centroid) == 0) {
      writeLines(paste("No contacts found for file:", path_aus), "errors.log", append = TRUE)
      return(list(ok = FALSE, filename = basename(path_aus), error = "No residues within cutoff distance"))
    }
    
    antigen_length <- length(unique(centroid_df$resno[centroid_df$chain == "A"])) # Differentiate between peptide-binding Ag and protein-binding Ag
    interface_radius <- round(max(vet_dist_centroid), 2)
    paratope_aa_number <- length(BS_HL)
    epitope_aa_number <- length(BS_A)
    interface_gyration_radius <- round(sqrt(mean(rowSums((df_centroids_BS[, c("x", "y", "z")] - center_BS)^2))), 2)
    
    # Add descriptors to output
    mean_NetDes_df$antigen_length <- antigen_length
    mean_NetDes_df$interface_radius <- interface_radius
    mean_NetDes_df$paratope_aa_number <- paratope_aa_number
    mean_NetDes_df$epitope_aa_number <- epitope_aa_number
    mean_NetDes_df$interface_gyration_radius <- interface_gyration_radius
    
    return(list(
      ok = TRUE,
      filename =pdbL[s],
      mean_NetDes = mean_NetDes_df
    ))
    
  }, error = function(e) {
    cat("Error in file:", pdbL[s], "\n", e$message, "\n")
    return(list(ok = FALSE, filename = pdbL[s], error = e$message))
  })
}

# Execute
with_progress({
  results_list <- future_map(1:length(pdbL), process_pdb_file, .progress = TRUE, .options = furrr_options(seed = TRUE))
})

valid_results <- results_list[map_lgl(results_list, "ok")]
failed_files <- map_chr(discard(results_list, "ok"), "filename")

# Build a unique dataframe
df_NetDes_general <- bind_rows(map(valid_results, "mean_NetDes"))

# Save as .csv
write.csv(df_NetDes_general, file = file.path(files_dir, "output_descrittori.csv"), row.names = FALSE)
