### THIS SCRIPT ANALYZES ANTIGEN–ANTIBODY INTERACTIONS FROM MULTIPLE PDB FILES ###
### IT CALCULATES RESIDUE–RESIDUE CONTACT MATRICES AND AMINO ACID COUNTS WITHIN CONCENTRIC "RINGS" AROUND THE INTERFACE CENTER ###
### RESULTS ARE SAVED AS CSV AND RDS FILES FOR EACH RING ###

### IMPORTANT CONVENTION FOR CONTACT MATRICES ###
# - In the **asymmetric contact matrices**:
#   * ROWS (dimension 1) correspond to AMINO ACIDS from the ANTIBODY (chains H/L)
#   * COLUMNS (dimension 2) correspond to AMINO ACIDS from the ANTIGEN (chain A)
#
# - Example:
#   If ALA from the antibody contacts ARG from the antigen,
#   the cell [ALA, ARG] (row ALA, column ARG) is incremented.
#
# - If the same amino acid appears in both chains (e.g. ALA Ab – ALA Ag),
#   the corresponding diagonal cell [ALA, ALA] is incremented
#   → this is still a valid cross-chain contact, not a self-contact within one chain.
#
# - In the **symmetric contact matrices**:
#   * Each contact is added twice, once in [aa1, aa2] and once in [aa2, aa1]
#   * Diagonal entries are incremented by +2 for each contact
#   * This makes the matrix symmetric by construction.
#
# => Always remember: in the asymmetric matrices,
#    ROWS = ANTIBODY (H/L), COLUMNS = ANTIGEN (A).

### Libraries ###
library(bio3d)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(latex2exp)
library(readr)

### Setting paths, directories and global variables ###
pdb_dir <- "/Users/lorenzosisti/Downloads/database_settembre_renamed/"
results_dir <- "/Users/lorenzosisti/Downloads/stats_int/"  
dir.create(results_dir, showWarnings = FALSE)

DistCutoff <- 8.5  # Two residues are defined as "in contact" if the centroids of their side chain are at less than 8.5 Angstroms
Nsteps <- 20 # Number of concentric "rings" around the interface

pdbL <- list.files(pdb_dir)

# Amino acids ordered according to the Kyte-Doolittle hydrophobicity scale
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE")

r_max <- 20
r_int_df <- data.frame(PDB = character(), r_int = numeric(), stringsAsFactors = FALSE)


### MAIN: GENERATING CONTACT MATRICES AND COMPUTING RESIDUE COUNTS ###

# Initialize cumulative data structures
residue_counts_interface_ab_sum <- vector("list", Nsteps) # residue counts for antibody chains (H/L)
residue_counts_interface_ligand_df_sum <- vector("list", Nsteps) # residue counts for antigen chain (A)

for (step in 1:Nsteps) {
  residue_counts_interface_ab_sum[[step]] <- setNames(rep(0, length(amino_acids)), amino_acids)
  residue_counts_interface_ligand_df_sum[[step]] <- setNames(rep(0, length(amino_acids)), amino_acids)
}


#Function to process each PDB file in the pdbL list
process_pdb_file <- function(s) {
  tryCatch({
    pdb_file <- pdbL[s]
    cat("Elaborating:", pdb_file, "\n")
    
    # Retrieve and read the PDB file
    path_aus <- file.path(pdb_dir, pdb_file)
    pdb_aus <- read.pdb(path_aus)
    file_name <- pdb_file
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
    # assegna nomi per matching (CRUCIALE per le estrazioni successive)
    if (nrow(df_centroids_BS) > 0) {
      names(vet_dist_centroid) <- rownames(df_centroids_BS)
    } else {
      vet_dist_centroid <- numeric(0)
    }
    
    # calcola r_int in modo robusto (NA se non ci sono centroids)
    if (length(vet_dist_centroid) > 0) {
      r_int <- max(vet_dist_centroid, na.rm = TRUE)
    } else {
      r_int <- NA_real_
      cat("Warning: r_int NA for", pdb_file, "- no interface centroids found\n", file = "errors.log", append = TRUE)
    }
    borders <- seq(0, r_max, length.out = Nsteps + 1)
    
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
    local_ab_counts <- vector("list", Nsteps)
    local_ag_counts <- vector("list", Nsteps)
    
    for (step in 1:Nsteps) {
      local_ab_counts[[step]] <- setNames(rep(0, length(amino_acids)), amino_acids)
      local_ag_counts[[step]] <- setNames(rep(0, length(amino_acids)), amino_acids)
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
      ab_counts = local_ab_counts,
      ag_counts = local_ag_counts,
      r_int = r_int,
      filename = pdb_file
      
    ))
    
    
    
  }, error = function(e) {
    cat("Error in", pdbL[s], ":", e$message, "\n", file = "errors.log", append = TRUE)
  })
}

# Execute
for (s in seq_along(pdbL)) {
  result <- process_pdb_file(s)
  if (is.null(result)) next  # Skip file if it returned NULL
  
  r_int_df <- rbind(r_int_df, data.frame(PDB = result$filename, r_int = result$r_int, stringsAsFactors = FALSE))
  
  for (step in 1:Nsteps) {
    residue_counts_interface_ab_sum[[step]] <- residue_counts_interface_ab_sum[[step]] + result$ab_counts[[step]]
    residue_counts_interface_ligand_df_sum[[step]] <- residue_counts_interface_ligand_df_sum[[step]] + result$ag_counts[[step]]
  }
}

# Save results 

write.csv(r_int_df, file.path(results_dir, "r_int_per_PDB.csv"), row.names = FALSE)

datasets <- list(
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

mean_r_int <- mean(r_int_df$r_int, na.rm = TRUE)

### CHECK ON THE CONTACT NUMBERS IN EACH RADIAL RING ###

### Load saved data ###
# leggi i conteggi pre-salvati
residue_counts_interface_ab <- readRDS(file.path(results_dir, "residue_counts_interface_ab.rds"))
residue_counts_interface_ligando <- readRDS(file.path(results_dir, "residue_counts_interface_ligand_df.rds"))

# prepara contenitori
combined_counts <- vector("list", Nsteps)   # somma ab + ag per ring
combined_freq   <- vector("list", Nsteps)   # frequenze normalizzate per ring

for (step in 1:Nsteps) {
  # Somma semplice dei conteggi (assumi che le due liste abbiano vettori nominati con gli stessi aminoacidi)
  counts_ab <- residue_counts_interface_ab[[step]]
  counts_ag <- residue_counts_interface_ligando[[step]]
  
  # Se per qualche motivo uno dei vettori non è presente, crealo a zero (robustezza)
  if (is.null(counts_ab)) counts_ab <- setNames(rep(0, length(amino_acids)), amino_acids)
  if (is.null(counts_ag)) counts_ag <- setNames(rep(0, length(amino_acids)), amino_acids)
  
  # Somma elemento per elemento (mantiene i nomi)
  combined <- counts_ab + counts_ag
  combined_counts[[step]] <- combined
  
  # Frequenze normalizzate: dividi per la somma totale del ring (se somma 0 -> frequenze tutte 0)
  total_combined <- sum(combined)
  if (total_combined > 0) {
    freq <- combined / total_combined
  } else {
    freq <- setNames(rep(0, length(combined)), names(combined))
  }
  combined_freq[[step]] <- freq
  
  # salva su disco (opzionale): counts e freq per ring
  write.csv(as.data.frame(combined), file = file.path(results_dir, paste0("combined_counts_step", step, ".csv")), row.names = TRUE)
  write.csv(as.data.frame(freq),     file = file.path(results_dir, paste0("combined_freq_step", step, ".csv")),   row.names = TRUE)
}

charged <- c("ARG", "LYS", "ASP", "GLU", "HIS")
polar <- c("ASN", "GLN", "SER", "THR")
hydrophobics <- c("PRO", "GLY", "ALA", "MET", "CYS", "LEU", "VAL", "ILE")
aromatics <- c("TYR", "TRP", "PHE")

# contenitori per risultati aggregati
group_names <- c("charged","polar","hydrophobic","aromatic","other")
group_freq_list <- vector("list", Nsteps)         # lista di vettori con 4/5 frequenze
group_freq_matrix <- matrix(0, nrow = length(group_names), ncol = Nsteps,
                            dimnames = list(group_names, paste0("step", 1:Nsteps)))

for(step in 1:Nsteps) {
  freq_vec <- combined_freq[[step]]
  # assicurati che il vettore abbia tutti i 20 aminoacidi nello stesso ordine
  freq_vec <- freq_vec[amino_acids]   # se manca, produce NA -> sostituisci con 0
  freq_vec[is.na(freq_vec)] <- 0
  
  # somma per gruppo
  s_charged <- sum(freq_vec[intersect(names(freq_vec), charged)])
  s_polar   <- sum(freq_vec[intersect(names(freq_vec), polar)])
  s_hydrophobic <- sum(freq_vec[intersect(names(freq_vec), hydrophobics)])
  s_aromatic <- sum(freq_vec[intersect(names(freq_vec), aromatics)])
  s_other <- 1 - (s_charged + s_polar + s_hydrophobic + s_aromatic)  # garantisce che i totali sommino a 1
  
  # se preferisci non forzare la somma a 1 (ma usare la somma reale), usa:
  # s_other <- sum(freq_vec[setdiff(names(freq_vec), all_groups)])
  
  # popola output
  group_freq_list[[step]] <- c(charged = s_charged,
                               polar = s_polar,
                               hydrophobic = s_hydrophobic,
                               aromatic = s_aromatic,
                               other = s_other)
  group_freq_matrix[, step] <- group_freq_list[[step]]
}

# controlli rapidi
rowSums(group_freq_matrix)   # dovrebbero essere tutti = 1 (o molto prossimi)
print(group_freq_matrix)

# salva su disco (opzionale)
write.csv(group_freq_matrix, file = file.path(results_dir, "group_freq_matrix.csv"), row.names = TRUE)

freq_df <- as.data.frame(group_freq_matrix)
freq_df <- freq_df[rownames(freq_df) != "other", ]

# Trasforma la matrice in formato long per ggplot2
# Trasforma i rownames in una colonna "Group"
freq_df$Group <- rownames(freq_df)

# Melt in formato long
freq_long <- melt(freq_df, id.vars = "Group", variable.name = "Step", value.name = "Frequency")

x_labels <- paste0(0:(length(freq_sub_t$Step)-1), "-", 1:length(freq_sub_t$Step))
my_labels <- c("charged_delta" = "Charged residues", "hydrophobic_delta" = "Hydrophobic residues")

# Grouped bar plot
ggplot(freq_long, aes(x = Step, y = Frequency, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  theme_minimal() +
  labs(title = "Frequenze amminoacidi per gruppo e step",
       x = "Distance from centroid", y = "Frequency") +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Colori e etichette dei gruppi
my_colors <- c("charged" = "gold1",
               "polar" = "mediumvioletred",
               "hydrophobic" = "dodgerblue2",
               "aromatic" = "grey40")

# Nuovi nomi per la legenda (in LaTeX)
my_labels <- c("charged" = TeX("Charged residues"),
               "polar" = TeX("Polar residues"),
               "hydrophobic" = TeX("Hydrophobic residues"),
               "aromatic" = TeX("Aromatic residues"))

# Etichette sull'asse x (0-1, 1-2, ... 19-20)
x_labels <- paste0(0:(ncol(freq_df)-2), "-", 1:(ncol(freq_df)-1))

# Trasforma in formato long e aggiungi Step come fattore
freq_long$Step <- factor(freq_long$Step, levels = colnames(freq_df)[1:(ncol(freq_df)-1)])

# Plot
ggplot(freq_long, aes(x = Step, y = Frequency, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  theme_minimal() +
  labs(#title = TeX("Frequenze amminoacidi per gruppo e step"),
       x = TeX("Distance from centroid (Å)"),
       y = TeX("Frequency"),
       fill = NULL) +
  scale_fill_manual(values = my_colors, labels = my_labels) +
  scale_x_discrete(labels = x_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        panel.grid.major.x = element_blank(),   # rimuove linee verticali
        panel.grid.minor = element_blank(),     # rimuove linee secondarie
        panel.grid.major.y = element_line(color = "grey80")) 
# Plot con barre più sottili
ggplot(freq_long, aes(x = Step, y = Frequency, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.75) +  # width < 1 rende le barre più sottili
  theme_minimal() +
  labs(#title = TeX("Frequenze amminoacidi per gruppo e step"),
    x = TeX("Distance from centroid (Å)"),
    y = TeX("Frequency"),
    fill = NULL) +
  scale_fill_manual(values = my_colors, labels = my_labels) +
  scale_x_discrete(labels = x_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        panel.grid.major.x = element_blank(),   
        panel.grid.minor = element_blank(),     
        panel.grid.major.y = element_line(color = "grey80"))




# Prendi solo le categorie di interesse
freq_sub <- freq_df[c("charged", "hydrophobic"), ]

# Trasponi per avere gli step sulle righe
freq_sub_t <- t(freq_sub)
freq_sub_t <- as.data.frame(freq_sub_t)
freq_sub_t$Step <- rownames(freq_sub_t)  # aggiungi colonna Step
# Rimuovi l'ultima riga
freq_sub_t <- freq_sub_t[-nrow(freq_sub_t), ]
# Calcola la differenza tra idrofobici e carichi
freq_sub_t$hydrophobic <- as.numeric(freq_sub_t$hydrophobic)
freq_sub_t$charged     <- as.numeric(freq_sub_t$charged)

freq_sub_t$charged_delta <- freq_sub_t$charged - min(freq_sub_t$charged)
freq_sub_t$hydrophobic_delta <- freq_sub_t$hydrophobic - min(freq_sub_t$hydrophobic)

# Trasforma in formato long
freq_long <- freq_sub_t %>%
  select(Step, charged_delta, hydrophobic_delta) %>%
  pivot_longer(cols = c(charged_delta, hydrophobic_delta),
               names_to = "Group", values_to = "Value")

# Ordina Step
freq_long$Step <- factor(freq_long$Step, levels = freq_sub_t$Step)

# Plot con curve smussate
ggplot(freq_long, aes(x = as.numeric(Step), y = Value, color = Group)) +
  #geom_line(size = 1.2) +
  #geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, span = 0.3) +  # curve smussate tipo densità
  theme_minimal() +
  labs(title = "Delta Frequenze: Hydrophobic vs Charged",
       x = "Step (distance from centroid)",
       y = "Delta Frequency",
       color = "Group") +
  scale_x_continuous(breaks = 1:length(freq_sub_t$Step), labels = freq_sub_t$Step) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Simula piccoli jitter attorno ai valori per creare "distribuzioni"
set.seed(123)
freq_long$Value_jitter <- freq_long$Value + rnorm(nrow(freq_long), mean = 0, sd = 0.005)

ggplot(freq_long, aes(x = Step, y = Value_jitter, color = Group, group = Group)) +
  geom_density(stat = "identity", position = "identity") # ma attenzione: con pochi punti il risultato non sarà affidabile

# Colori personalizzati
my_colors <- c("charged_delta" = "gold1", "hydrophobic_delta" = "dodgerblue2")

# Plot con area sottesa
ggplot(freq_long, aes(x = as.numeric(Step), y = Value, fill = Group)) +
  geom_area(position = "identity", alpha = 0.3, color = NA) +  # riempimento semitrasparente
  geom_line(method = "loess", se = FALSE, span = 0.3) +  # curve smussate tipo densità
  theme_minimal() +
  labs(title = "Delta Frequenze: Hydrophobic vs Charged",
       x = "Step (distance from centroid)",
       y = "Delta Frequency",
       color = "Group",
       fill = "Group") +
  scale_x_continuous(breaks = 1:length(freq_sub_t$Step), labels = freq_sub_t$Step) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Plot con area sottesa
ggplot(freq_long, aes(x = as.numeric(Step), y = Value, fill = Group)) +
  geom_area(position = "identity", alpha = 0.3, color = NA) +  # riempimento semitrasparente
  geom_line(aes(color = Group), size = 0.7) +                  # linee sopra il riempimento
  theme_minimal() +
  labs(title = "Delta Frequenze: Hydrophobic vs Charged",
       x = "Step (distance from centroid)",
       y = "Delta Frequency",
       color = "Group",
       fill = "Group") +
  scale_x_continuous(breaks = 1:length(freq_sub_t$Step), labels = freq_sub_t$Step) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Plot con area sottesa e legenda in alto
ggplot(freq_long, aes(x = as.numeric(Step), y = Value, fill = Group)) +
  geom_area(position = "identity", alpha = 0.3, color = NA) +  # riempimento semitrasparente
  geom_line(aes(color = Group), size = 0.7) +                  # linee sopra il riempimento
  theme_minimal() +
  labs(#title = TeX("Delta Frequenze: Hydrophobic vs Charged"),
       x = TeX("Distance from centroid (Å)"),
       y = TeX("$\\Delta$ Frequency"),
       color = "Group",
       fill = "Group") +
  scale_x_continuous(breaks = 1:length(freq_sub_t$Step), labels = x_labels) +
  scale_color_manual(values = my_colors, labels = my_labels) +
  scale_fill_manual(values = my_colors, labels = my_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        text = element_text(size = 12),
        panel.grid = element_blank())

residue_counts_paratope <- readRDS("/Users/lorenzosisti/Downloads/matrice_singola_settembre/residue_counts_interface_ab_sum.rds")
residue_counts_epitope <- readRDS("/Users/lorenzosisti/Downloads/matrice_singola_settembre/residue_counts_interface_ligando_sum.rds")

paratope_freq <- residue_counts_paratope / sum(residue_counts_paratope)
epitope_freq <- residue_counts_epitope / sum(residue_counts_epitope)

delta_freq <- paratope_freq - epitope_freq

delta_df <- data.frame(
  AminoAcid = names(delta_freq),
  DeltaFreq = as.numeric(delta_freq)
)

# Ordina per arricchimento/deplezione
delta_df <- delta_df[order(delta_df$DeltaFreq, decreasing = TRUE), ]

# Plot
ggplot(delta_df, aes(x = reorder(AminoAcid, DeltaFreq), y = DeltaFreq, fill = DeltaFreq > 0)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(x = "Amino Acid", y = "ΔΔ Frequency (Paratope - Epitope)",
       title = "Differenza composizione paratope vs epitope") +
  scale_fill_manual(values = c("TRUE" = "dodgerblue2", "FALSE" = "firebrick3"),
                    labels = c("Epitope enriched", "Paratope enriched"),
                    guide = "legend")

ggplot(delta_df, aes(x = reorder(AminoAcid, DeltaFreq), y = DeltaFreq, fill = DeltaFreq > 0)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    x = TeX("Amino Acid"),
    y = TeX("$\\Delta$ Frequency (Paratope - Epitope)"),
    title = TeX("Composition difference at the Paratope-Epitope interface"),
    fill = NULL
  ) +
  scale_fill_manual(
    values = c("TRUE" = "dodgerblue2", "FALSE" = "gold1"),
    labels = c("Epitope-enriched amino acids", "Paratope-enriched amino acids")
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.x = element_blank(),   # niente linee verticali
    panel.grid.minor = element_blank(),     # niente linee secondarie
    panel.grid.major.y = element_line(color = "grey80"), # solo linee orizzontali
    plot.title = element_text(hjust = 0.5)  # titolo centrato
  )


# Carica i dati
descriptors <- read.csv("/Users/lorenzosisti/Downloads/interface_descriptors_nuovo/merged_descriptors.csv")

# Seleziona le colonne
selected_cols <- descriptors[, c("epitope_aa_number", "paratope_aa_number")]

# Calcola la differenza
selected_cols$difference <- selected_cols$paratope_aa_number - selected_cols$epitope_aa_number
selected_cols$sum <- selected_cols$paratope_aa_number + selected_cols$epitope_aa_number

mean_aa_numb_int <- mean(selected_cols$sum)

# Visualizza la differenza
head(selected_cols)

# Statistiche di base sulla differenza
summary(selected_cols$difference)

# Istogramma della differenza
ggplot(selected_cols, aes(x = difference)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Distribuzione delle differenze tra paratope e epitope",
    x = "Differenza (Paratope - Epitope)",
    y = "Frequenza"
  )


# Calcola le categorie
selected_cols <- selected_cols %>%
  mutate(category = case_when(
    difference >= -3 & difference <= 3 ~ "Comparable patches",
    difference > 3 ~ "Paratope patch >> Epitope patch",
    difference < -3 ~ "Epitope patch >> Paratope patch"
  ))

# Conta le occorrenze per categoria
pie_data <- selected_cols %>%
  group_by(category) %>%
  summarise(count = n()) %>%
  mutate(perc = count / sum(count) * 100)

# Grafico a torta
ggplot(pie_data, aes(x = "", y = perc, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  geom_text(aes(label = paste0(round(perc, 1), "%")),
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("Comparable patches" = "gold1",
                               "Paratope >> Epitope" = "dodgerblue",
                               "Epitope >> Paratope" = "mediumvioletred")) +
  labs(fill = "Category", title = "Distribuzione delle differenze tra paratope e epitope")

library(ggplot2)
library(latex2exp)

# Grafico a torta con LaTeX
ggplot(pie_data, aes(x = "", y = perc, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  geom_text(aes(label = paste0(round(perc, 1), "%")),
            position = position_stack(vjust = 0.5), size = 5) +
  scale_fill_manual(values = c("Comparable patches" = "gold1",
                               "Paratope patch >> Epitope patch" = "dodgerblue",
                               "Epitope patch >> Paratope patch" = "mediumvioletred")) +
  labs(fill = NULL, 
       title = NULL)

