#########################################################################################
# THIS SCRIPT ANALYZES ANTIGEN-ANTIBODY INTERACTIONS FROM MULTIPLE PDB FILES.
# FOR EACH PDB IT:
# 1) COMPUTES SIDE-CHAIN CENTROIDS (VIA renumber_ab_chains() IN functions.R)
# 2) IDENTIFIES THE INTERFACE RESIDUES (ANTIBODY H/L vs ANTIGEN = ALL OTHER CHAINS)
# 3) ASSIGNS EACH INTERFACE RESIDUE TO ONE OF Nsteps CONCENTRIC "RINGS" AROUND THE
#    INTERFACE CENTER (CENTROID OF ALL INTERFACE RESIDUES), UP TO r_max ANGSTROM
# 4) COUNTS THE AMINO ACIDS PER RING, SEPARATELY FOR ANTIBODY AND ANTIGEN
# RESULTS (SUMMED OVER ALL PDBs) ARE SAVED AS CSV AND RDS FILES FOR EACH RING.
#
# CONVENTIONS:
# - Antibody = chains H and L; antigen = ALL the other chains.
# - Rings are [a, b): residues farther than r_max from the interface center
#   are not counted (their number is reported in n_beyond_rmax).
# - r_int (interface radius) is computed on ALL interface residues,
#   including those beyond r_max.
#########################################################################################

### Libraries ###
pacman::p_load(bio3d, dplyr, future, furrr, purrr, progressr, reshape2, latex2exp)
source("/Users/lorenzosisti/TiNDER/functions.R")

### Set up parallelization ###
plan(multisession, workers = parallel::detectCores() - 1)
handlers(global = TRUE)
handlers("rstudio")

### Setting paths, directories and global variables ###
pdb_dir <- "/Users/lorenzosisti/Downloads/database_settembre_renamed/"
results_dir <- "/Users/lorenzosisti/TiNDER/data/interface_statistics_fig_2/"
dir.create(results_dir, showWarnings = FALSE)

DistCutoff <- 8.5  # Two residues are "in contact" if the centroids of their side chains are at less than 8.5 Angstrom
Nsteps <- 20       # Number of concentric "rings" around the interface
r_max <- 20        # Outer radius (Angstrom) of the last ring

# Amino acids ordered according to the Kyte-Doolittle hydrophobicity scale
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE")

pdb_paths <- list.files(pdb_dir, pattern = "\\.pdb$", full.names = TRUE)
borders <- seq(0, r_max, length.out = Nsteps + 1)


### Process a single PDB file ###
### N.B. This function relies on the renumber_ab_chains() function (see functions.R)
process_pdb_file <- function(pdb_path) {
  file_name <- basename(pdb_path)
  tryCatch({
    pdb_aus <- read.pdb(pdb_path)
    chain_HL <- c("H", "L")
    
    renumbered_df <- renumber_ab_chains(pdb_aus, pdb_path, log_file = "errors.log")
    if (!renumbered_df$ok) {
      return(list(ok = FALSE, filename = file_name, path = pdb_path, error = renumbered_df$error))
    }
    
    # Centroid of each residue
    centroids_df <- as.data.frame(
      renumbered_df$df_coord_renumbered %>%
        group_by(chain, resno, resid) %>%
        summarise(x = mean(x), y = mean(y), z = mean(z), .groups = 'drop')
    )
    centroids_df$residue_key <- paste(centroids_df$resid, centroids_df$resno, centroids_df$chain, sep = "_")
    rownames(centroids_df) <- centroids_df$residue_key
    
    # Antibody-antigen distances (drop = FALSE avoids collapsing to a vector)
    DistMat <- as.matrix(dist(centroids_df[, c("x", "y", "z")]))
    condHL <- centroids_df$chain %in% chain_HL
    condA  <- !condHL
    Inter_DistMat_Bin <- DistMat[condA, condHL, drop = FALSE] <= DistCutoff
    
    BS_A_names  <- rownames(Inter_DistMat_Bin)[rowSums(Inter_DistMat_Bin) > 0]
    BS_HL_names <- colnames(Inter_DistMat_Bin)[colSums(Inter_DistMat_Bin) > 0]
    if (length(BS_A_names) == 0 || length(BS_HL_names) == 0) {
      return(list(ok = FALSE, filename = file_name, path = pdb_path, error = "no contacts"))
    }
    
    # Interface center and distance of each interface residue from it
    df_BS <- centroids_df[c(BS_A_names, BS_HL_names), ]
    center_BS <- colMeans(df_BS[, c("x", "y", "z")])
    df_BS$dist_centroid <- sqrt((df_BS$x - center_BS["x"])^2 +
                                  (df_BS$y - center_BS["y"])^2 +
                                  (df_BS$z - center_BS["z"])^2)
    r_int <- max(df_BS$dist_centroid)
    
    # Ring assignment: 1..Nsteps; beyond r_max -> Nsteps + 1 (not counted)
    df_BS$ring <- findInterval(df_BS$dist_centroid, borders)
    n_beyond_rmax <- sum(df_BS$ring > Nsteps)
    
    # Amino acid counts per ring (antibody = H/L, antigen = all other chains)
    count_ring <- function(df, step) {
      setNames(as.numeric(table(factor(df$resid[df$ring == step], levels = amino_acids))),
               amino_acids)
    }
    df_ab <- df_BS[df_BS$chain %in% chain_HL, ]
    df_ag <- df_BS[!df_BS$chain %in% chain_HL, ]
    
    list(
      ok = TRUE,
      filename = file_name,
      r_int = r_int,
      n_beyond_rmax = n_beyond_rmax,
      ab_counts = lapply(1:Nsteps, function(step) count_ring(df_ab, step)),
      ag_counts = lapply(1:Nsteps, function(step) count_ring(df_ag, step))
    )
  }, error = function(e) {
    list(ok = FALSE, filename = file_name, path = pdb_path, error = e$message)
  })
}


### MAIN: run in parallel ###
with_progress({
  results_list <- future_map(
    pdb_paths,
    process_pdb_file,
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
  )
})

### Valid results and failed files ###
valid_results <- keep(results_list, ~ .x$ok)
failed_files  <- map_chr(discard(results_list, ~ .x$ok), "filename")
cat("File falliti:", length(failed_files), "\n")
if (length(failed_files) > 0) writeLines(failed_files, file.path(results_dir, "failed_files.txt"))

### Interface radius per PDB ###
r_int_df <- map_dfr(valid_results, ~ data.frame(PDB = .x$filename, r_int = .x$r_int,
                                                n_beyond_rmax = .x$n_beyond_rmax))
mean_r_int <- mean(r_int_df$r_int, na.rm = TRUE)

### Residue counts per ring, summed over all PDBs ###
residue_counts_interface_ab_sum <- lapply(1:Nsteps, function(step)
  Reduce(`+`, map(valid_results, ~ .x$ab_counts[[step]])))
residue_counts_interface_ligand_sum <- lapply(1:Nsteps, function(step)
  Reduce(`+`, map(valid_results, ~ .x$ag_counts[[step]])))


### Save results ###
write.csv(r_int_df, file.path(results_dir, "r_int_per_PDB.csv"), row.names = FALSE)

datasets <- list(
  residue_counts_interface_ab = residue_counts_interface_ab_sum,
  residue_counts_interface_ligand = residue_counts_interface_ligand_sum
)

for (name in names(datasets)) {
  for (step in 1:Nsteps) {
    fname <- file.path(results_dir, paste0(name, "_step", step, ".csv"))
    write.csv(datasets[[name]][[step]], file = fname, row.names = TRUE)
  }
  saveRDS(datasets[[name]], file = file.path(results_dir, paste0(name, ".rds")))
}
#########################################################################################
# 1) COMPOSITION OF THE INTERFACE IN EACH RADIAL RING (amino acid groups, charged vs hydrophobic)
# 2) COMPOSITION DIFFERENCE BETWEEN PARATOPE AND EPITOPE (per amino acid)
# 3) SIZE OF PARATOPE vs EPITOPE (two overlaid PDFs)
#########################################################################################

### Libraries ###
pacman::p_load(dplyr, tidyr, ggplot2, latex2exp)

### Parameters (the same of the ring script: redundant if run in the same session) ###
results_dir <- "/Users/lorenzosisti/TiNDER/data/interface_statistics_fig_2/"
Nsteps <- 20   # Number of concentric rings (1 Angstrom each, since r_max = 20)
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE")

# Common theme for all the plots (same base_size -> same font size everywhere)
tema_pub <- theme_minimal() +
  theme(legend.position = "top",
        panel.grid.minor = element_blank())

# Ring labels on the x axis (0-1, 1-2, ... 19-20)
x_labels <- paste0(0:(Nsteps - 1), "-", 1:Nsteps)


### Load saved data (lists of Nsteps named vectors, one per ring) ###
residue_counts_interface_ab     <- readRDS(file.path(results_dir, "residue_counts_interface_ab.rds"))
residue_counts_interface_ligand <- readRDS(file.path(results_dir, "residue_counts_interface_ligand.rds"))


################################################################################
### 1) COMPOSITION IN EACH RING
################################################################################

# Sum of antibody + antigen counts and normalized frequencies, per ring
combined_counts <- vector("list", Nsteps)
combined_freq   <- vector("list", Nsteps)

for (step in 1:Nsteps) {
  combined <- residue_counts_interface_ab[[step]] + residue_counts_interface_ligand[[step]]
  combined_counts[[step]] <- combined
  
  # If the ring is empty -> frequencies all 0
  freq <- if (sum(combined) > 0) combined / sum(combined) else combined
  combined_freq[[step]] <- freq
  
  write.csv(as.data.frame(combined), file = file.path(results_dir, paste0("combined_counts_step", step, ".csv")), row.names = TRUE)
  write.csv(as.data.frame(freq),     file = file.path(results_dir, paste0("combined_freq_step", step, ".csv")),   row.names = TRUE)
}

# The four groups cover all the 20 amino acids (no "other" group needed)
groups <- list(
  charged     = c("ARG", "LYS", "ASP", "GLU", "HIS"),
  polar       = c("ASN", "GLN", "SER", "THR"),
  hydrophobic = c("PRO", "GLY", "ALA", "MET", "CYS", "LEU", "VAL", "ILE"),
  aromatic    = c("TYR", "TRP", "PHE")
)

# Rows = groups, columns = rings
group_freq_matrix <- sapply(combined_freq, function(freq_vec) sapply(groups, function(g) sum(freq_vec[g])))
colnames(group_freq_matrix) <- paste0("step", 1:Nsteps)

# Quick check: each ring (column) should sum to 1
colSums(group_freq_matrix)
print(group_freq_matrix)
write.csv(group_freq_matrix, file = file.path(results_dir, "group_freq_matrix.csv"), row.names = TRUE)

my_colors <- c("charged" = "gold1",
               "polar" = "mediumvioletred",
               "hydrophobic" = "dodgerblue2",
               "aromatic" = "grey40")

labels_groups <- c("charged" = TeX("Charged residues"),
                   "polar" = TeX("Polar residues"),
                   "hydrophobic" = TeX("Hydrophobic residues"),
                   "aromatic" = TeX("Aromatic residues"))

### Grouped bar plot: frequency of each amino acid group per ring ###
freq_long <- group_freq_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Group") %>%
  pivot_longer(-Group, names_to = "Step", values_to = "Frequency") %>%
  mutate(Step = factor(Step, levels = paste0("step", 1:Nsteps)),
         Group = factor(Group, levels = names(groups)))

p_group_bars <- ggplot(freq_long, aes(x = Step, y = Frequency, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.75) +
  labs(x = TeX("Distance from centroid (Å)"),
       y = TeX("Frequency"),
       fill = NULL) +
  scale_fill_manual(values = my_colors, labels = labels_groups) +
  scale_x_discrete(labels = x_labels) +
  tema_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank())

print(p_group_bars)

### Heatmap: amino acid groups (rows) x rings (columns), plasma palette ###
labels_groups_plain <- c("charged" = "Charged", "polar" = "Polar",
                         "hydrophobic" = "Hydrophobic", "aromatic" = "Aromatic")

p_group_heatmap <- ggplot(freq_long, aes(x = Step,
                                         y = factor(Group, levels = rev(names(groups))),
                                         fill = Frequency)) +
  geom_tile(color = "white", linewidth = 0.5) +
  labs(x = TeX("Distance from interface centroid (Å)"),
       y = NULL,
       fill = "Frequency") +
  scale_x_discrete(labels = x_labels) +
  scale_y_discrete(labels = labels_groups_plain) +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  tema_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(hjust = 0),
        panel.grid = element_blank())

print(p_group_heatmap)

### Charged vs hydrophobic: frequency relative to the minimum value across rings ###
delta_long <- as.data.frame(t(group_freq_matrix[c("charged", "hydrophobic"), ])) %>%
  mutate(Step = 1:Nsteps,
         charged_delta = charged - min(charged),
         hydrophobic_delta = hydrophobic - min(hydrophobic)) %>%
  pivot_longer(c(charged_delta, hydrophobic_delta), names_to = "Group", values_to = "Value")

colors_delta <- c("charged_delta" = "purple", "hydrophobic_delta" = "#CC4678FF")
labels_delta <- c("charged_delta" = TeX("Charged residues"),
                  "hydrophobic_delta" = TeX("Hydrophobic residues"))

p_delta_rings <- ggplot(delta_long, aes(x = Step, y = Value, fill = Group)) +
  geom_area(position = "identity", alpha = 0.3, color = NA) +
  geom_line(aes(color = Group), linewidth = 0.7) +
  labs(x = TeX("Distance from interface centroid (Å)"),
       y = TeX("$\\Delta$ Frequency"),
       color = NULL, fill = NULL) +
  scale_x_continuous(breaks = 1:Nsteps, labels = x_labels) +
  scale_color_manual(values = colors_delta, labels = labels_delta) +
  scale_fill_manual(values = colors_delta, labels = labels_delta) +
  tema_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

print(p_delta_rings)

################################################################################
### 2) COMPOSITION DIFFERENCE: PARATOPE vs EPITOPE (per amino acid)
################################################################################

# Counts summed over all rings (same dataset and same r_max as the ring analysis)
paratope_counts <- Reduce(`+`, residue_counts_interface_ab)
epitope_counts  <- Reduce(`+`, residue_counts_interface_ligand)

paratope_freq <- paratope_counts / sum(paratope_counts)
epitope_freq  <- epitope_counts / sum(epitope_counts)

delta_df <- data.frame(
  AminoAcid = names(paratope_freq),
  DeltaFreq = as.numeric(paratope_freq - epitope_freq)
) %>%
  mutate(Enriched = ifelse(DeltaFreq > 0, "Paratope", "Epitope"))

# One-letter code of the amino acids (shown in the central band of the plot)
aa_one_letter <- c(ARG = "R", LYS = "K", ASN = "N", ASP = "D", GLN = "Q", GLU = "E",
                   HIS = "H", PRO = "P", TYR = "Y", TRP = "W", SER = "S", THR = "T",
                   GLY = "G", ALA = "A", MET = "M", CYS = "C", PHE = "F", LEU = "L",
                   VAL = "V", ILE = "I")

# Half-height of the central band that hosts the letters (in frequency units).
# The bars start at +-gap and the y axis labels are shifted accordingly,
# so the axis still shows the true values of DeltaFreq.
gap <- 0.012

# Left -> right: from the most paratope-enriched to the most epitope-enriched amino acid
delta_plot_df <- delta_df %>%
  arrange(desc(DeltaFreq)) %>%
  mutate(x = row_number(),
         letter = unname(aa_one_letter[AminoAcid]),
         ymin = sign(DeltaFreq) * gap,
         ymax = DeltaFreq + sign(DeltaFreq) * gap)

y_breaks_true <- pretty(delta_df$DeltaFreq, n = 6)
y_breaks_true <- y_breaks_true[y_breaks_true != 0]
y_breaks_pos  <- y_breaks_true + sign(y_breaks_true) * gap

p_delta_aa <- ggplot(delta_plot_df) +
  geom_hline(yintercept = c(-gap, gap), color = "grey80", linewidth = 0.3) +
  geom_rect(aes(xmin = x - 0.4, xmax = x + 0.4, ymin = ymin, ymax = ymax,
                fill = Enriched, color = Enriched),
            alpha = 0.4, linewidth = 0.8) +
  geom_text(aes(x = x, y = 0, label = letter), fontface = "bold", size = 4.5) +
  labs(
    x = NULL,
    y = TeX("$\\Delta$ Frequency (Paratope - Epitope)"),
    title = TeX("Composition difference at the Paratope-Epitope interface"),
    fill = NULL, color = NULL
  ) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = y_breaks_pos, labels = y_breaks_true) +
  scale_fill_manual(
    values = c("Paratope" = "#CC4678FF", "Epitope" = "purple"),
    labels = c("Paratope" = "Paratope-enriched amino acids", "Epitope" = "Epitope-enriched amino acids")
  ) +
  scale_color_manual(
    values = c("Paratope" = "#CC4678FF", "Epitope" = "purple"),
    labels = c("Paratope" = "Paratope-enriched amino acids", "Epitope" = "Epitope-enriched amino acids")
  ) +
  tema_pub +
  theme(panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

print(p_delta_aa)


################################################################################
### 3) SIZE OF PARATOPE vs EPITOPE: TWO OVERLAID PDFs
################################################################################

descriptors <- read.csv("/Users/lorenzosisti/Downloads/interface_descriptors_nuovo/merged_descriptors.csv")

# Mean number of interface amino acids (paratope + epitope)
mean_aa_numb_int <- mean(descriptors$paratope_aa_number + descriptors$epitope_aa_number)

size_long <- descriptors %>%
  select(paratope_aa_number, epitope_aa_number) %>%
  pivot_longer(everything(), names_to = "Patch", values_to = "n_aa")

colors_patch <- c("paratope_aa_number" = "#CC4678FF", "epitope_aa_number" = "purple")
labels_patch <- c("paratope_aa_number" = "Paratope", "epitope_aa_number" = "Epitope")

p_size_pdf <- ggplot(size_long, aes(x = n_aa, fill = Patch, color = Patch)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  labs(
    x = "Number of amino acids in the interface patch",
    y = "Probability Density Function (PDF)",
    fill = NULL, color = NULL
  ) +
  scale_fill_manual(values = colors_patch, labels = labels_patch) +
  scale_color_manual(values = colors_patch, labels = labels_patch) +
  tema_pub

print(p_size_pdf)
