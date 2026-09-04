# --- 1. CARICAMENTO LIBRERIE ---
library(dplyr)
library(ggplot2)
library(ggExtra)
library(patchwork)
library(pROC)
library(tidyr)

# --- 2. IMPOSTAZIONI GLOBALI ---
# !!! MODIFICA QUESTO PATH !!!
output_dir <- "/Users/lorenzosisti/Downloads/DockQ_HDOCK_decoy_vs_all_09_08"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
plot_width <- 8
plot_height <- 6
plot_dpi <- 300
theme_custom <- theme_minimal() +
  theme(plot.title = element_blank())

path_dockq_af3 = "/Users/lorenzosisti/Downloads/DockQ_results_AF3_12_06.csv"
path_pot_af3_whole = "/Users/lorenzosisti/Downloads/potenziali_statistici_30_06_50_pose/potenziali_whole_per_posa.csv"
path_pot_af3_strat = "/Users/lorenzosisti/Downloads/potenziali_statistici_30_06_50_pose/potenziali_gr_per_posa.csv"
path_pot_af3_cdr = "/Users/lorenzosisti/Downloads/potenziali_statistici_30_06_50_pose/punteggi_cdr_per_posa.csv"

#path_dockq_af3 = "/Users/lorenzosisti/Downloads/DockQ_HDOCK/DockQ_results_HDOCK.csv"
#path_pot_af3_whole = "/Users/lorenzosisti/Downloads/potenziali_statistici_30_06_hdock/potenziali_whole_per_posa.csv"
#path_pot_af3_strat = "/Users/lorenzosisti/Downloads/potenziali_statistici_30_06_hdock/potenziali_gr_per_posa.csv"
#path_pot_af3_cdr = "/Users/lorenzosisti/Downloads/potenziali_statistici_30_06_hdock/punteggi_cdr_per_posa.csv"

dockq_scores <- read.csv(path_dockq_af3)
whole_potential_scores <- read.csv(path_pot_af3_whole)
stratified_potential_scores <- read.csv(path_pot_af3_strat)
cdr_potentials_scores <- read.csv(path_pot_af3_cdr)

# --- 3. CALCOLO DELLE 6 ROC (ogni oggetto ha ora un nome distinto) ---

## --- WHOLE ---
whole_merged_df <- dockq_scores %>%
  left_join(whole_potential_scores, by = c("Model" = "pdb"))
df_roc_whole <- whole_merged_df
df_roc_whole$true_class <- ifelse(df_roc_whole$DockQ <= 0.24, 1, 0)

roc_whole_sym <- roc(df_roc_whole$true_class, df_roc_whole$mean_sym)
print(roc_whole_sym)
roc_whole_asym <- roc(df_roc_whole$true_class, df_roc_whole$mean_asym)
print(roc_whole_asym)

## --- STRATIFIED ---
strat_merged_df <- dockq_scores %>%
  left_join(stratified_potential_scores, by = c("Model" = "pdb"))
df_roc_strat <- strat_merged_df
df_roc_strat$true_class <- ifelse(df_roc_strat$DockQ <= 0.24, 1, 0)

roc_strat_sym <- roc(df_roc_strat$true_class, df_roc_strat$mean_sym)
print(roc_strat_sym)
roc_strat_asym <- roc(df_roc_strat$true_class, df_roc_strat$mean_asym)
print(roc_strat_asym)

## --- CDR ---
cdr_merged_df <- dockq_scores %>%
  left_join(cdr_potentials_scores, by = c("Model" = "pdb_filename"))
df_roc_cdr <- cdr_merged_df
df_roc_cdr$true_class <- ifelse(df_roc_cdr$DockQ <= 0.24, 1, 0)

roc_cdr_sym <- roc(df_roc_cdr$true_class, df_roc_cdr$mean_sym)
print(roc_cdr_sym)
roc_cdr_asym <- roc(df_roc_cdr$true_class, df_roc_cdr$mean_asym)
print(roc_cdr_asym)

# --- 4. PLOT COMBINATO DELLE 6 ROC NELLO STESSO GRAFICO ---

roc_list <- list(
  "Whole - Sym"      = roc_whole_sym,
  "Whole - Asym"     = roc_whole_asym,
  "Stratified - Sym" = roc_strat_sym,
  "Stratified - Asym"= roc_strat_asym,
  "CDR - Sym"        = roc_cdr_sym,
  "CDR - Asym"       = roc_cdr_asym
)

# Aggiungo l'AUC di ciascuna curva direttamente nell'etichetta di legenda
auc_values <- sapply(roc_list, function(r) sprintf("%.3f", as.numeric(auc(r))))
names(roc_list) <- paste0(names(roc_list), " (AUC = ", auc_values, ")")

p_roc_combined <- ggroc(roc_list, legacy.axes = TRUE, linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  labs(
    x = "1 - Specificità (False Positive Rate)",
    y = "Sensitività (True Positive Rate)",
    color = "Curva ROC"
  ) +
  theme_custom +
  theme(legend.position = "right")

print(p_roc_combined)

ggsave(
  filename = file.path(output_dir, "AF3_ROC_combined_6curve_mean.png"),
  plot = p_roc_combined,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi
)
