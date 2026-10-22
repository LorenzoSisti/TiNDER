#########################################################################################
# CORRELATION BETWEEN THE UPPER AND THE LOWER TRIANGLE OF THE ASYMMETRIC POTENTIAL MATRIX
# 1) Builds the 20x20 asymmetric potential matrix (rows = antibody, columns = antigen),
#    with the amino acids ordered according to the Kyte-Doolittle scale
# 2) For each pair of amino acids (i, j), i != j, compares
#      upper triangle: V(Ab = i, Ag = j)   vs   lower triangle: V(Ab = j, Ag = i)
#    (190 pairs; the diagonal is excluded: it is the same element in both triangles)
# 3) Scatterplot of the two triangles, each point colored by the hydropathy
#    complementarity Hr of the amino acid pair (Grassmann et al., JCIM 2025)
#
# N.B. In the Hr scale low values = strong complementarity, values close to 1 = weak.
#      With the default plasma direction, purple = low Hr = high hydropathy complementarity.
#########################################################################################

pacman::p_load(readr, dplyr, tidyr, ggplot2)

# Amino acids ordered according to the Kyte-Doolittle hydrophobicity scale
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE")

hydropathy_complementarity_matrix <- read_csv("/Users/lorenzosisti/TiNDER/data/hydropathy_complementarity/Hr_complementarity_matrix.csv")
statistical_potential_matrix <- read_csv("/Users/lorenzosisti/TiNDER/data/pmf_whole_cdr_sippl_18_10/sc_centroid/whole_int_asym_potential.csv")

whole_interface_statistical_potential_matrix <- statistical_potential_matrix[statistical_potential_matrix$part == "all", ]


### 1) 20x20 asymmetric potential matrix (rows = resid_ab, columns = resid_ag) ###
potential_wide <- pivot_wider(whole_interface_statistical_potential_matrix,
                              id_cols = resid_ab, names_from = resid_ag, values_from = potential)
potential_matrix <- as.matrix(potential_wide[, -1])
rownames(potential_matrix) <- potential_wide$resid_ab
potential_matrix <- potential_matrix[amino_acids, amino_acids]

# Hydropathy complementarity matrix (first column of the csv = amino acid names)
Hr_matrix <- as.matrix(hydropathy_complementarity_matrix[, -1])
rownames(Hr_matrix) <- hydropathy_complementarity_matrix[[1]]
Hr_matrix <- Hr_matrix[amino_acids, amino_acids]


### 2) Upper vs lower triangle, pair by pair ###
upper_idx <- which(upper.tri(potential_matrix), arr.ind = TRUE)   # (row i, col j) with i < j

pairs_df <- data.frame(
  aa_i  = amino_acids[upper_idx[, "row"]],
  aa_j  = amino_acids[upper_idx[, "col"]],
  upper = potential_matrix[upper_idx],         # V(Ab = i, Ag = j)
  lower = potential_matrix[upper_idx[, 2:1]],  # V(Ab = j, Ag = i)
  # Hr of the unordered pair (mean of the two entries, robust even if Hr were not perfectly symmetric)
  Hr    = (Hr_matrix[upper_idx] + Hr_matrix[upper_idx[, 2:1]]) / 2
)

r_triangles <- cor(pairs_df$upper, pairs_df$lower, method = "pearson")
cat("Pearson correlation between upper and lower triangle:", round(r_triangles, 3),
    "(n =", nrow(pairs_df), "pairs)\n")


### 3) Scatterplot colored by hydropathy complementarity ###
p_triangles <- ggplot(pairs_df, aes(x = upper, y = lower, color = Hr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5, alpha = 0.9) +
  annotate("text", 
           x = mean(range(pairs_df$upper)), # Punto medio orizzontale
           y = Inf,                        # Bordo superiore
           hjust = 0.5,                    # Centrato orizzontalmente
           vjust = 1,                    # Più in alto
           label = sprintf("Pearson r = %.2f", r_triangles)) +
  scale_color_viridis_c(option = "plasma") +   # purple = low Hr = high hydropathy complementarity
  labs(x = "Upper triangle potential (Ab residue i - Ag residue j)",
       y = "Lower triangle potential (Ab residue j - Ag residue i)",
       color = expression(H[r])) +
  coord_equal() +
  theme_minimal() +
  theme(legend.title = element_text(hjust = 0.2)) # 0.5 = centro, aumenta a 0.7 o 1 per spostarlo più a destra

print(p_triangles)
