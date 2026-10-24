### FIG_3_C SCRIPT

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(viridis)

output_dir <- "/Users/lorenzosisti/TiNDER/data/Fig_3_C"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

potentials_cor <- read_csv("/Users/lorenzosisti/TiNDER/data/pmf_whole_cdr_sippl_24_10/sc_centroid/cor_asym_potential.csv")

pair_cor <- potentials_cor %>%
  pivot_longer(
    cols = -part,
    names_to = "part2",
    values_to = "correlation"
  ) %>%
  filter(
    part != part2,
    correlation <= 0.4,
    match(part, c("all", "h1", "h2", "h3", "l1", "l2", "l3")) <
      match(part2, c("all", "h1", "h2", "h3", "l1", "l2", "l3"))
  ) %>%
  select(part, part2, correlation)

potentials_matrices <- read_csv("/Users/lorenzosisti/TiNDER/data/pmf_whole_cdr_sippl_24_10/sc_centroid/whole_int_asym_potential.csv")

residuals_all <- pair_cor %>%
  rename(cor_matrix = correlation) %>%
  pmap_dfr(function(part, part2, cor_matrix){
    
    mat1 <- potentials_matrices %>%
      filter(.data$part == part) %>%
      select(resid_ag, resid_ab, potential) %>%
      rename(potential_x = potential)
    
    mat2 <- potentials_matrices %>%
      filter(.data$part == part2) %>%
      select(resid_ag, resid_ab, potential) %>%
      rename(potential_y = potential)
    
    df <- inner_join(mat1, mat2,
                     by = c("resid_ag","resid_ab"))
    
    fit <- lm(potential_y ~ potential_x, data = df)
    
    df <- df %>%
      mutate(
        part = part,
        part2 = part2,
        fitted = fitted(fit),
        residual = residuals(fit),
        distance = abs(residual)
      )
    
    r <- cor(df$potential_x, df$potential_y)
    
    p <- ggplot(
      df,
      aes(x = potential_x, y = potential_y, color = distance)
    ) +
      geom_smooth(
        aes(x = potential_x, y = potential_y),
        method = "lm",
        formula = y ~ x,
        se = FALSE,
        color = "grey40",
        linewidth = 0.8,
        inherit.aes = FALSE
      ) +
      geom_point(size = 2.5, alpha = 0.9) +
      annotate(
        "text",
        x = mean(range(df$potential_x)),
        y = Inf,
        hjust = 0.5,
        vjust = 1.3,
        label = sprintf("Pearson r = %.2f", r)
      ) +
      scale_color_viridis_c(option = "viridis") +
      labs(
        title = paste(part, "vs", part2),
        x = paste0(part, " potential"),
        y = paste0(part2, " potential"),
        color = "Absolute residual"
      ) +
      coord_equal() +
      theme_minimal() +
      theme(legend.title = element_text(hjust = 0.5))
    
    ggsave(
      filename = file.path(output_dir,
                           paste0("scatter_", part, "_vs_", part2, ".pdf")),
      plot = p,
      width = 5.5,
      height = 5.5
    )
    
    df
  })

residuals_dist <- residuals_all %>%
  filter(is.finite(distance))

p_residuals <- ggplot(residuals_dist, aes(x = distance)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 60,
    fill = "grey75",
    color = "white",
    alpha = 0.8
  ) +
  geom_density(
    color = viridis::viridis(1),
    linewidth = 1
  ) +
  labs(
    title = "Distribution of absolute residuals",
    subtitle = "All selected pairs of potential matrices",
    x = "Absolute residual",
    y = "Density"
  ) +
  theme_minimal(base_size = 12)

print(p_residuals)

ggsave(
  filename = file.path(
    output_dir,
    "absolute_residuals_distribution.pdf"
  ),
  plot = p_residuals,
  width = 7,
  height = 5
)


residuals_dist %>%
  summarise(
    n = n(),
    median = median(distance),
    q90 = quantile(distance, 0.90),
    q95 = quantile(distance, 0.95),
    q99 = quantile(distance, 0.99),
    q995 = quantile(distance, 0.995),
    q999 = quantile(distance, 0.999),
    maximum = max(distance)
  )

#threshold <- 1.48
threshold <- 0.807

pairs_above_threshold <- residuals_dist %>%
  filter(distance > threshold) %>%
  select(
    part,
    part2,
    resid_ag,
    resid_ab,
    potential_x,
    potential_y,
    fitted,
    residual,
    distance
  ) %>%
  arrange(desc(distance))

# Controllo del risultato
dim(pairs_above_threshold)
head(pairs_above_threshold)

# Esportazione in CSV
write.csv(
  pairs_above_threshold,
  file = file.path(output_dir, "pairs_above_q999_residual.csv"),
  row.names = FALSE
)

plot_data <- pairs_above_threshold %>%
  mutate(
    amino_acid_pair = paste(resid_ag, resid_ab, sep = " - "),
    matrix_pair = paste(part, part2, sep = " vs "),
    label = paste0(amino_acid_pair, "\n(", matrix_pair, ")")
  ) %>%
  # Conta quante volte compare ogni coppia di residui
  add_count(resid_ag, resid_ab, name = "n_observations") %>%
  arrange(desc(distance)) %>%
  mutate(
    observation = row_number(),
    label = factor(label, levels = unique(label))
  )

p_outliers <- ggplot(
  plot_data,
  aes(x = observation, y = distance, fill = n_observations)
) +
  geom_col(width = 0.8) +
  scale_fill_viridis_c(
    option = "plasma",
    name = "Number of\nobservations"
  ) +
  scale_x_continuous(
    breaks = plot_data$observation,
    labels = plot_data$label
  ) +
  labs(
    title = "Residue pairs with large absolute residuals",
    subtitle = paste0(
      "Absolute residual > ", threshold,
      " | Sorted in decreasing order"
    ),
    x = "Amino acid pair (matrix comparison)",
    y = "Absolute residual from linear fit"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 7
    ),
    panel.grid.major.x = element_blank()
  )

print(p_outliers)

# Esporta in PDF
ggsave(
  filename = file.path(output_dir, "barplot_large_residuals.pdf"),
  plot = p_outliers,
  width = max(10, nrow(plot_data) * 0.25),
  height = 7,
  limitsize = FALSE
)
















# ============================================================
# FIG_3_C: ANALISI DEGLI ABSOLUTE RESIDUALS
# PER CLASSI CHIMICO-FISICHE DI CONTATTO
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

# ------------------------------------------------------------
# 1. DEFINIZIONE DELLE CLASSI AMMINOACIDICHE
# ------------------------------------------------------------

groups <- list(
  PC = c("ARG", "LYS", "HIS"), # positively charged
  NC = c("ASP", "GLU"), # negatively charged
  P = c("ASN", "GLN", "SER", "THR"), # polar
  H = c("PRO", "GLY", "ALA", "MET", "CYS", "LEU", "VAL", "ILE"), # hydrophobic
  A = c("TYR", "TRP", "PHE") # aromatic
)

class_levels <- names(groups)

# Funzione per assegnare la classe chimica.
# Accetta codici come ARG, ARG45, 45_ARG o A:ARG:45.

classify_aa <- function(x) {
  vapply(as.character(x), function(z) {
    
    z <- toupper(z)
    
    hits <- names(groups)[
      vapply(groups, function(aa) {
        any(vapply(
          aa,
          function(a) grepl(a, z, fixed = TRUE),
          logical(1)
        ))
      }, logical(1))
    ]
    
    if (length(hits) == 1) {
      hits
    } else {
      NA_character_
    }
    
  }, character(1))
}

# ------------------------------------------------------------
# 2. ASSEGNAZIONE DELLE CLASSI AI DUE RESIDUI
# ------------------------------------------------------------

residuals_classified <- residuals_dist %>%
  mutate(
    class_ag = classify_aa(resid_ag),
    class_ab = classify_aa(resid_ab)
  )

# Controlla eventuali residui non classificati
unclassified <- residuals_classified %>%
  filter(is.na(class_ag) | is.na(class_ab))

cat("Osservazioni totali:", nrow(residuals_classified), "\n")
cat("Osservazioni non classificabili:", nrow(unclassified), "\n")

# Mantiene soltanto le osservazioni con entrambi i residui classificati
residuals_classified <- residuals_classified %>%
  filter(!is.na(class_ag), !is.na(class_ab))

# ------------------------------------------------------------
# 3. CREAZIONE DELLE 15 CLASSI DI CONTATTO
#    A-B equivale a B-A
# ------------------------------------------------------------

residuals_classified <- residuals_classified %>%
  mutate(
    index_ag = match(class_ag, class_levels),
    index_ab = match(class_ab, class_levels),
    
    class_1 = if_else(
      index_ag <= index_ab, class_ag, class_ab
    ),
    
    class_2 = if_else(
      index_ag <= index_ab, class_ab, class_ag
    ),
    
    contact_class = paste(class_1, class_2, sep = " - ")
  )

contact_levels <- unlist(lapply(
  seq_along(class_levels),
  function(i) {
    paste(class_levels[i:length(class_levels)],
          collapse = "")
  }
))

# Ordine esplicito delle 15 combinazioni
contact_levels <- unlist(lapply(
  seq_along(class_levels),
  function(i) {
    vapply(
      i:length(class_levels),
      function(j) {
        paste(class_levels[i], class_levels[j], sep = " - ")
      },
      character(1)
    )
  }
))

residuals_classified <- residuals_classified %>%
  mutate(
    contact_class = factor(
      contact_class,
      levels = contact_levels
    )
  )

# ------------------------------------------------------------
# 4. STATISTICHE PER CIASCUNA CLASSE DI CONTATTO
# ------------------------------------------------------------

class_summary <- residuals_classified %>%
  group_by(class_1, class_2, contact_class) %>%
  summarise(
    n_total = n(),
    
    median_distance = median(distance),
    mean_distance = mean(distance),
    q90 = quantile(distance, 0.90),
    q95 = quantile(distance, 0.95),
    
    n_outliers = sum(distance > threshold),
    pct_outliers = 100 * mean(distance > threshold),
    
    median_outlier = if (any(distance > threshold)) {
      median(distance[distance > threshold])
    } else {
      NA_real_
    },
    
    .groups = "drop"
  ) %>%
  arrange(desc(pct_outliers), desc(median_distance))

print(class_summary)

write.csv(
  class_summary,
  file.path(output_dir, "residuals_by_contact_class.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. FIGURA 1: DISTRIBUZIONE DI TUTTI I RESIDUI
# ------------------------------------------------------------

p_all_residuals <- ggplot(
  residuals_classified,
  aes(x = contact_class, y = distance, fill = contact_class)
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    alpha = 0.8,
    color = "grey35"
  ) +
  geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    color = "grey20"
  ) +
  scale_fill_viridis_d(option = "plasma", guide = "none") +
  labs(
   # title = "Absolute residuals by chemical contact class",
   # subtitle = "All finite observations",
    x = "Chemical contact class",
    y = "Distance from linear regression fit"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      #angle = 60,
      hjust = 0.5,
      vjust = 1,
      size = 8
    ),
    panel.grid.major.x = element_blank()
  )

print(p_all_residuals)

ggsave(
  file.path(output_dir, "violin_all_residuals_by_class.pdf"),
  p_all_residuals,
  width = 12,
  height = 7
)

# ------------------------------------------------------------
# 6. FIGURA 2: HEATMAP DELLA PERCENTUALE DI OUTLIER
# ------------------------------------------------------------

# Statistiche per la heatmap: tutte le osservazioni classificate
heatmap_summary <- residuals_classified %>%
  group_by(class_1, class_2) %>%
  summarise(
    n_total = n(),
    n_outliers = sum(distance > threshold),
    pct_outliers = 100 * mean(distance > threshold),
    .groups = "drop"
  )

# Espande la tabella in una matrice simmetrica 5 x 5
heatmap_data <- bind_rows(
  heatmap_summary,
  heatmap_summary %>%
    transmute(
      class_1 = class_2,
      class_2 = class_1,
      n_total,
      n_outliers,
      pct_outliers
    )
) %>%
  distinct(class_1, class_2, .keep_all = TRUE) %>%
  complete(
    class_1 = class_levels,
    class_2 = class_levels,
    fill = list(
      n_total = 0L,
      n_outliers = 0L,
      pct_outliers = NA_real_
    )
  ) %>%
  mutate(
    class_1 = factor(class_1, levels = class_levels),
    class_2 = factor(class_2, levels = class_levels)
  )

p_heatmap <- ggplot(
  heatmap_data,
  aes(x = class_1, y = class_2, fill = pct_outliers)
) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(
    aes(
      label = ifelse(
        n_total > 0,
        sprintf(
          "%.1f%%\n%d/%d",
          pct_outliers, n_outliers, n_total
        ),
        ""
      )
    ),
    size = 3.1
  ) +
  scale_fill_viridis_c(
    option = "plasma",
    name = "Outliers (%)",
    na.value = "white"
  ) +
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(limits = rev(class_levels)) +
  coord_fixed() +
  labs(
    title = "Outlier frequency by chemical contact class",
    subtitle = paste0(
      "Absolute residual > ", threshold,
      " | Cell label: percentage and n outliers / n total"
    ),
    x = "First chemical class",
    y = "Second chemical class"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 0, hjust = 1,
    )
  )

print(p_heatmap)

ggsave(
  file.path(output_dir, "heatmap_outlier_frequency.pdf"),
  p_heatmap,
  width = 9,
  height = 8
)


# Barplot delle frequenze di outlier per classe di contatto

plot_frequency <- heatmap_summary %>%
  mutate(
    contact_class = paste(class_1, class_2, sep = " - "),
    label = sprintf("%.1f%%", pct_outliers)
  ) %>%
  filter(n_total > 0) %>%
  arrange(desc(pct_outliers)) %>%
  mutate(
    contact_class = factor(
      contact_class,
      levels = unique(contact_class)
    )
  )

p_frequency <- ggplot(
  plot_frequency,
  aes(
    x = contact_class,
    y = pct_outliers
  )
) +
  geom_col(
    fill = "#CC4678FF",
    width = 0.8
  ) +
  geom_text(
    aes(label = label),
    vjust = -0.4,
    size = 3.2
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    # title = "Outlier frequency by chemical contact class",
    # subtitle = paste0(
    #   "Absolute residual > ", threshold,
    #   " | Ordered by decreasing frequency"
    # ),
    x = "Chemical contact class",
    y = "Outliers (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      color = "black",
      #angle = 45,
      hjust = 0.5,
      vjust = 1
    ),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

print(p_frequency)

ggsave(
  file.path(output_dir, "barplot_outlier_frequency.pdf"),
  p_frequency,
  width = 12,
  height = 7,
  limitsize = FALSE
)
