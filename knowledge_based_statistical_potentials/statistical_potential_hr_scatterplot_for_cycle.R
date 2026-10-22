pacman::p_load(readr, dplyr, tidyr, ggplot2)

# Amminoacidi ordinati secondo la scala Kyte-Doolittle
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE")

# Caricamento dataset
hydropathy_complementarity_matrix <- read_csv("/Users/lorenzosisti/TiNDER/data/hydropathy_complementarity/Hr_complementarity_matrix.csv")
#statistical_potential_matrix <- read_csv("/Users/lorenzosisti/TiNDER/data/pmf_whole_cdr_sippl_18_10/sc_centroid/whole_int_asym_potential.csv") %>%
  #mutate(part = toupper(part))
statistical_potential_matrix <- read_csv("/Users/lorenzosisti/TiNDER/data/pmf_layer_sippl_18_10/ring_asym_potential.csv") %>%
  rename(part = ring_pair) %>%
  mutate(part = toupper(part))

output_dir <- "/Users/lorenzosisti/TiNDER/data/scatterplots_gr"
# Crea la cartella se non esiste già
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 1) PRE-ELABORAZIONE (eseguita una sola volta fuori dal ciclo)
Hr_matrix <- as.matrix(hydropathy_complementarity_matrix[, -1])
rownames(Hr_matrix) <- hydropathy_complementarity_matrix[[1]]
Hr_matrix <- Hr_matrix[amino_acids, amino_acids]

# Indici del triangolo superiore (20x20 costante)
upper_idx <- which(upper.tri(Hr_matrix), arr.ind = TRUE)

# Calcolo del valore Hr per ogni coppia (invariante rispetto a 'part')
Hr_values <- (Hr_matrix[upper_idx] + Hr_matrix[upper_idx[, 2:1]]) / 2

# Lista di tutte le 'part' presenti nel file
unique_parts <- unique(toupper(statistical_potential_matrix$part))

# Range globale del potenziale (su tutte le part, prima del ciclo)
global_range <- range(statistical_potential_matrix$potential, na.rm = TRUE)
pad <- diff(global_range) * 0.05
global_lims <- c(global_range[1] - pad, global_range[2] + pad)

# 2) CICLO FOR SULLE DIVERSE 'PART'
for (p in unique_parts) {
  
  # Filtraggio e pivoting per la singola part
  potential_sub <- statistical_potential_matrix %>% 
    filter(part == p)
  
  potential_wide <- pivot_wider(
    potential_sub,
    id_cols = resid_ab,
    names_from = resid_ag,
    values_from = potential
  )
  
  potential_matrix <- as.matrix(potential_wide[, -1])
  rownames(potential_matrix) <- potential_wide$resid_ab
  potential_matrix <- potential_matrix[amino_acids, amino_acids]
  
  # Costruzione data frame triangolo superiore vs inferiore
  pairs_df <- data.frame(
    aa_i  = amino_acids[upper_idx[, "row"]],
    aa_j  = amino_acids[upper_idx[, "col"]],
    upper = potential_matrix[upper_idx],
    lower = potential_matrix[upper_idx[, 2:1]],
    Hr    = Hr_values
  )
  
  r_triangles <- cor(pairs_df$upper, pairs_df$lower, method = "pearson", use = "complete.obs")
  cat(sprintf("Part: %-12s | Pearson r = %.3f\n", p, r_triangles))
  
  # Generazione dello scatterplot
  p_triangles <- ggplot(pairs_df, aes(x = upper, y = lower, color = Hr)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 2.5, alpha = 0.9) +
    annotate(
      "text", 
      x = mean(range(pairs_df$upper, na.rm = TRUE)), 
      y = Inf, 
      hjust = 0.5, 
      vjust = 1.2, 
      label = sprintf("Pearson r = %.2f", r_triangles)
    ) +
    scale_color_viridis_c(option = "plasma") +
    coord_fixed(ratio = 1, xlim = global_lims, ylim = global_lims, expand = FALSE) +  # <-- qui
    labs(
      title = paste("Contacting layers:", p), # Here, put "CDR" when working w/ CDRs and not layers
      x = "Upper triangle potential (Ab residue i - Ag residue j)",
      y = "Lower triangle potential (Ab residue j - Ag residue i)",
      color = expression(H[r])
    ) +
    theme_minimal() +
    theme(legend.title = element_text(hjust = 0.2))
  
  # Salvataggio del grafico in PNG (dimensioni 7x5 o 5x7 pollici a seconda dell'orientamento desiderato)
  output_filename <- paste0("scatterplot_potential_", p, ".png")
  ggsave(
    filename = output_filename,
    path = output_dir,        
    plot = p_triangles,
    width = 7,
    height = 5,
    dpi = 300
  )
}
