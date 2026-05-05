# --- 1. CARICAMENTO LIBRERIE ---
library(dplyr)
library(ggplot2)
library(ggExtra)
library(patchwork)
library(pROC)
library(tidyr) 

# --- 2. IMPOSTAZIONI GLOBALI ---

# !!! MODIFICA QUESTO PATH !!!
output_dir <- "/Users/lorenzosisti/Downloads/potenziali_statistici_AF3_prova_strati_opt_maggio/ROC" 

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

plot_width <- 8
plot_height <- 6
plot_dpi <- 300

theme_custom <- theme_minimal() +
  theme(plot.title = element_blank())

#path_dockq_af3 = "/Users/lorenzosisti/Downloads/DockQ_HDOCK/DockQ_results_AF3.csv"
#path_pot_af3_whole = "/Users/lorenzosisti/Downloads/potenziali_statistici_AF3_opt/potenziali_per_modello_AF3.csv"
#path_pot_af3_strat = "/Users/lorenzosisti/Downloads/potenziali_statistici_AF3_prova_strati_opt_aprile/summary_results_AF3_stratificati_sparse.csv"

path_dockq_af3 = "/Users/lorenzosisti/Downloads/DockQ_HDOCK/DockQ_results_HDOCK.csv"
path_pot_af3_whole = "/Users/lorenzosisti/Downloads/potenziali_statistici_HDOCK_opt/potenziali_per_modello_HDOCK.csv"
path_pot_af3_strat = "/Users/lorenzosisti/Downloads/potenziali_statistici_HDOCK_prova_strati_opt/summary_results_HDOCK_stratificati_sparse.csv"


#data("aSAH")
#roc1 <- roc(aSAH$outcome, aSAH$s100b)
#print(roc1)

dockq_scores <- read.csv(path_dockq_af3)
whole_potential_scores <- read.csv(path_pot_af3_whole)
stratified_potential_scores <- read.csv(path_pot_af3_strat)

whole_merged_df <- dockq_scores %>%
  left_join(whole_potential_scores, by = c("Model" = "pdb"))

df_roc_whole <- whole_merged_df[whole_merged_df$DockQ <= 0.24 | whole_merged_df$DockQ >= 0.81, ]
df_roc_whole$true_class <- ifelse(df_roc_whole$DockQ >= 0.81, 1, 0)

whole_sym_roc <- roc(df_roc_whole$true_class, df_roc_whole$sum_sym)
print(whole_sym_roc)

whole_asym_roc <- roc(df_roc_whole$true_class, df_roc_whole$sum_asym)
print(whole_asym_roc)

###

whole_merged_df <- dockq_scores %>%
  left_join(stratified_potential_scores, by = c("Model" = "pdb"))

df_roc_whole <- whole_merged_df[whole_merged_df$DockQ <= 0.24 | whole_merged_df$DockQ >= 0.81, ]
df_roc_whole$true_class <- ifelse(df_roc_whole$DockQ >= 0.81, 1, 0)

whole_sym_roc <- roc(df_roc_whole$true_class, df_roc_whole$sum_sym)
print(whole_sym_roc)

whole_asym_roc <- roc(df_roc_whole$true_class, df_roc_whole$sum_asym)
print(whole_asym_roc)
