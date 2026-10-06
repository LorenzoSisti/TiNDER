# --- 1. CARICAMENTO LIBRERIE ---
library(dplyr)
library(ggplot2)
library(ggExtra)
library(patchwork)
library(pROC)
library(tidyr)

# --- 2. IMPOSTAZIONI GLOBALI ---
# !!! MODIFICA QUESTO PATH !!!
output_dir <- "/Users/lorenzosisti/TiNDER/data/roc_auc_decoy_vs_all_sum_pmf"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
plot_width <- 8
plot_height <- 6
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

# --- ISTOGRAMMA DI fnat ---

p_fnat <- ggplot(dockq_scores, aes(x = fnat)) +
  geom_histogram(
    binwidth = 0.05,
    boundary = 0,          # allinea i bin a 0, 0.05, 0.10, ...
    closed   = "left",
    fill     = "dodgerblue2",
    color    = "white"
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    limits = c(-0.025, 1.025)
  ) +
  labs(
    x = "fnat (fraction of native contacts)",
    y = "Docking model counts"
  ) +
  theme_custom

print(p_fnat)

# --- 3. CALCOLO DELLE 6 ROC (ogni oggetto ha ora un nome distinto) ---

## --- WHOLE ---
whole_merged_df <- dockq_scores %>%
  left_join(whole_potential_scores, by = c("Model" = "pdb"))
df_roc_whole <- whole_merged_df
df_roc_whole$true_class <- ifelse(df_roc_whole$DockQ <= 0.24, 1, 0)

roc_whole_sym <- roc(df_roc_whole$true_class, df_roc_whole$sum_sym)
print(roc_whole_sym)
roc_whole_asym <- roc(df_roc_whole$true_class, df_roc_whole$sum_asym)
print(roc_whole_asym)

## --- STRATIFIED ---
strat_merged_df <- dockq_scores %>%
  left_join(stratified_potential_scores, by = c("Model" = "pdb"))
df_roc_strat <- strat_merged_df
df_roc_strat$true_class <- ifelse(df_roc_strat$DockQ <= 0.24, 1, 0)

roc_strat_sym <- roc(df_roc_strat$true_class, df_roc_strat$sum_sym)
print(roc_strat_sym)
roc_strat_asym <- roc(df_roc_strat$true_class, df_roc_strat$sum_asym)
print(roc_strat_asym)

## --- CDR ---
cdr_merged_df <- dockq_scores %>%
  left_join(cdr_potentials_scores, by = c("Model" = "pdb_filename"))
df_roc_cdr <- cdr_merged_df
df_roc_cdr$true_class <- ifelse(df_roc_cdr$DockQ <= 0.24, 1, 0)

roc_cdr_sym <- roc(df_roc_cdr$true_class, df_roc_cdr$score_global_sym)
print(roc_cdr_sym)
roc_cdr_asym <- roc(df_roc_cdr$true_class, df_roc_cdr$score_global_asym)
print(roc_cdr_asym)

# --- 4. PLOT COMBINATO DELLE 6 ROC NELLO STESSO GRAFICO ---

roc_list <- list(
  "Whole - Sym"       = roc_whole_sym,
  "Whole - Asym"      = roc_whole_asym,
  "Stratified - Sym"  = roc_strat_sym,
  "Stratified - Asym" = roc_strat_asym,
  "CDR - Sym"         = roc_cdr_sym,
  "CDR - Asym"        = roc_cdr_asym
)

# Palette personalizzata: ciano, gold, magenta (scuro = Sym, chiaro = Asym)
my_colors <- c(
  "Whole - Sym"       = "dodgerblue4",  # ciano scuro
  "Whole - Asym"      = "dodgerblue2",  # ciano chiaro
  "Stratified - Sym"  = "darkgoldenrod",  # gold scuro
  "Stratified - Asym" = "darkgoldenrod1",  # gold chiaro
  "CDR - Sym"         = "deeppink4",  # magenta scuro
  "CDR - Asym"        = "deeppink"   # magenta chiaro
)

# Etichette con AUC
auc_values <- sapply(roc_list, function(r) sprintf("%.3f", as.numeric(auc(r))))
new_labels <- paste0(names(roc_list), " (AUC = ", auc_values, ")")

# Rinomino lista e palette con le stesse etichette
names(roc_list)  <- new_labels
names(my_colors) <- new_labels

p_roc_combined <- ggroc(roc_list, legacy.axes = TRUE, linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = my_colors, breaks = new_labels) +
  labs(
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
    color = "ROC AUC values"
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
