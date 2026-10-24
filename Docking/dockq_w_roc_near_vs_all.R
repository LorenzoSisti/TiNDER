# --- 1. CARICAMENTO LIBRERIE ---
library(dplyr)
library(ggplot2)
library(ggExtra)
library(patchwork)
library(pROC)
library(tidyr)
library(viridisLite)

# --- 2. IMPOSTAZIONI GLOBALI ---
# !!! MODIFICA QUESTO PATH !!!
output_dir <- "/Users/lorenzosisti/TiNDER/data/roc_auc_near_vs_all_mean_pmf"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
plot_width <- 7
plot_height <- 5
plot_dpi <- 300
theme_custom <- theme_minimal() +
  theme(plot.title = element_blank())

path_dockq_af3 = "/Users/lorenzosisti/TiNDER/data/dockq_af3.csv"
path_pot_af3_whole = "/Users/lorenzosisti/TiNDER/data/score_per_af3_pose/whole_interface_score.csv"
path_pot_af3_strat = "/Users/lorenzosisti/TiNDER/data/score_per_af3_pose/layer_score.csv"
path_pot_af3_cdr = "/Users/lorenzosisti/TiNDER/data/score_per_af3_pose/cdr_score.csv"

dockq_scores <- read.csv(path_dockq_af3)
whole_potential_scores <- read.csv(path_pot_af3_whole)
stratified_potential_scores <- read.csv(path_pot_af3_strat)
cdr_potentials_scores <- read.csv(path_pot_af3_cdr)

# --- 3. CALCOLO DELLE 6 ROC (ogni oggetto ha ora un nome distinto) ---

## --- WHOLE ---
whole_merged_df <- dockq_scores %>%
  left_join(whole_potential_scores, by = c("Model" = "pdb"))
df_roc_whole <- whole_merged_df
df_roc_whole$true_class <- ifelse(df_roc_whole$DockQ <= 0.81, 1, 0)

roc_whole_sym <- roc(df_roc_whole$true_class, df_roc_whole$mean_sym)
print(roc_whole_sym)
roc_whole_asym <- roc(df_roc_whole$true_class, df_roc_whole$mean_asym)
print(roc_whole_asym)

## --- STRATIFIED ---
strat_merged_df <- dockq_scores %>%
  left_join(stratified_potential_scores, by = c("Model" = "pdb"))
df_roc_strat <- strat_merged_df
df_roc_strat$true_class <- ifelse(df_roc_strat$DockQ <= 0.81, 1, 0)

roc_strat_sym <- roc(df_roc_strat$true_class, df_roc_strat$mean_sym)
print(roc_strat_sym)
roc_strat_asym <- roc(df_roc_strat$true_class, df_roc_strat$mean_asym)
print(roc_strat_asym)

## --- CDR ---
cdr_merged_df <- dockq_scores %>%
  left_join(cdr_potentials_scores, by = c("Model" = "pdb_filename"))
df_roc_cdr <- cdr_merged_df
df_roc_cdr$true_class <- ifelse(df_roc_cdr$DockQ <= 0.81, 1, 0)

roc_cdr_sym <- roc(df_roc_cdr$true_class, df_roc_cdr$mean_sym)
print(roc_cdr_sym)
roc_cdr_asym <- roc(df_roc_cdr$true_class, df_roc_cdr$mean_asym)
print(roc_cdr_asym)

# --- 4. PLOT COMBINATO DELLE 6 ROC NELLO STESSO GRAFICO ---

roc_list <- list(
  "Whole-interface - Sym"       = roc_whole_sym,
  "Whole-interface - Asym"      = roc_whole_asym,
  "Layer-based - Sym"  = roc_strat_sym,
  "Layer-based - Asym" = roc_strat_asym,
  "CDR-based - Sym"         = roc_cdr_sym,
  "CDR-based - Asym"        = roc_cdr_asym
)

# Palette viridisLite: viola = AUC maggiore, giallo = AUC minore
auc_num   <- sapply(roc_list, function(r) as.numeric(auc(r)))
#my_colors <- viridisLite::plasma(6)[rank(-auc_num)]

# Palette plasma: esclude i toni più chiari iniziali
my_palette <- viridisLite::plasma(
  6,
  begin = 0,
  end = 0.8
)

my_colors <- my_palette[rank(-auc_num)]
# Etichette con AUC
#auc_values <- sprintf("%.3f", auc_num)
#new_labels <- paste0(names(roc_list), " (AUC = ", auc_values, ")")
new_labels <- names(roc_list)


# Rinomino lista e palette con le stesse etichette
names(roc_list)  <- new_labels
names(my_colors) <- new_labels

p_roc_combined <- ggroc(roc_list, legacy.axes = TRUE, linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = my_colors,
                     breaks = new_labels[order(auc_num, decreasing = TRUE)]) +
  labs(
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
    color = "ROC AUC"
  ) +
  theme_custom +
  theme(legend.position = "right")

print(p_roc_combined)

ggsave(
  filename = file.path(output_dir, "af3_roc_mean_pmf.png"),
  plot = p_roc_combined,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi
)

