# --- 1. CARICAMENTO LIBRERIE ---
library(dplyr)
library(pROC)
library(ggplot2)

# --- 2. PERCORSI FILE ---
path_dockq_af3 <- "/Users/lorenzosisti/Downloads/DockQ_HDOCK/DockQ_results_AF3.csv"
path_pot_af3_whole <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali_marzo/risultati_scoring_finale_AF3_WholeInt/MASTER_RESULTS_AF3_tidy.csv"
path_pot_af3_strat <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali_marzo/risultati_scoring_finale_AF3/MASTER_RESULTS_AF3_tidy.csv"

# --- 3. FUNZIONE DI CALCOLO ROC AUC RAGGRUPPATA ---
calcola_statistiche_auc <- function(dockq_path, pot_path) {
  
  # Caricamento dei file
  dockq_df <- read.csv(dockq_path)
  pot_df <- read.csv(pot_path)
  
  # Unione dei dataframe usando "Model" e "pdb"
  merged_df <- dockq_df %>%
    inner_join(pot_df, by = c("Model" = "pdb"))
  
  # 1. Filtro pose (rimuovo intermedie) e creo true_class binaria
  df_roc <- merged_df %>%
    filter(DockQ <= 0.24 | DockQ >= 0.81) %>%
    mutate(true_class = ifelse(DockQ >= 0.81, 1, 0))
  
  # 2. Calcolo della ROC AUC per ogni singola combinazione
  auc_per_group <- df_roc %>%
    group_by(seed, split_type, group_num, potential_type) %>%
    summarise(
      # Per calcolare la ROC ci servono almeno una posa positiva e una negativa nel gruppo
      n_classes = n_distinct(true_class),
      roc_auc = ifelse(n_classes == 2, 
                       as.numeric(roc(true_class, sum_potential, quiet = TRUE)$auc), 
                       NA),
      .groups = "drop"
    ) %>%
    # Rimuoviamo i gruppi dove l'AUC non è calcolabile (es. il gruppo aveva solo pose sbagliate)
    filter(!is.na(roc_auc)) 
  
  # 3. Media e Deviazione Standard per ogni split e tipo di potenziale (sym/asym)
  auc_summary <- auc_per_group %>%
    group_by(split_type, potential_type) %>%
    summarise(
      mean_auc = mean(roc_auc, na.rm = TRUE),
      sd_auc = sd(roc_auc, na.rm = TRUE),
      n_valid_groups = n(), # Quanti gruppi sono stati usati per questa media
      .groups = "drop"
    )
  
  return(list(
    dettaglio_gruppi = auc_per_group,
    riassunto_split = auc_summary
  ))
}

# --- 4. ESECUZIONE SUI DUE DATASET ---

# Esecuzione per Dataset Whole Interface
risultati_whole <- calcola_statistiche_auc(path_dockq_af3, path_pot_af3_whole)
cat("\n=== Risultati Whole Interface ===\n")
print(risultati_whole$riassunto_split)

# Esecuzione per Dataset Stratified
risultati_strat <- calcola_statistiche_auc(path_dockq_af3, path_pot_af3_strat)
cat("\n=== Risultati Stratified ===\n")
print(risultati_strat$riassunto_split)

# --- 5. SALVATAGGIO IN CSV ---
# !!! INSERISCI QUI IL PERCORSO DELLA TUA CARTELLA !!!
output_dir_saturazione <- "/Users/lorenzosisti/Downloads/saturazione_prestazione_potenziali_marzo/plot_saturazione_potenziali"
dir.create(output_dir_saturazione, showWarnings = FALSE, recursive = TRUE)

# Salvataggio dei riassunti in CSV
write.csv(risultati_whole$riassunto_split, 
          file.path(output_dir_saturazione, "ROC_AUC_Whole_summary.csv"), 
          row.names = FALSE)
write.csv(risultati_strat$riassunto_split, 
          file.path(output_dir_saturazione, "ROC_AUC_Stratified_summary.csv"), 
          row.names = FALSE)

# --- 6. PREPARAZIONE DATI PER IL PLOT DI SATURAZIONE ---
# Uniamo i due dataset aggiungendo una colonna identificativa
df_whole <- risultati_whole$riassunto_split %>% mutate(Dataset = "Whole Interface")
df_strat <- risultati_strat$riassunto_split %>% mutate(Dataset = "Stratified")
df_plot <- bind_rows(df_whole, df_strat)

# Estraiamo il numero di split dalla stringa e calcoliamo i PDB
df_plot <- df_plot %>%
  mutate(
    # Estrae solo le cifre dalla stringa "split_X"
    n_split = as.numeric(gsub("split_", "", split_type)),
    # Calcola il numero di PDB usati per l'allenamento
    num_pdbs = 2187 / n_split
  )

# --- 7. PLOT DELLA CURVA DI SATURAZIONE ---
# Creiamo il grafico con le curve per 'asym' e 'sym', separate per Dataset
p_saturazione <- ggplot(df_plot, aes(x = num_pdbs, y = mean_auc, color = potential_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  # Aggiungiamo le barre di errore per la deviazione standard
  geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc), 
                width = max(df_plot$num_pdbs) * 0.02, # Larghezza barre proporzionale all'asse x
                alpha = 0.6) +
  # Dividiamo il grafico in due pannelli: uno per Whole, uno per Stratified
  facet_wrap(~ Dataset) +
  theme_minimal() +
  scale_color_manual(values = c("sym" = "dodgerblue2", "asym" = "mediumvioletred")) +
  labs(
    x = "Numero stimato di PDB (2187 / n_split)",
    y = "Mean ROC AUC (\u00B1 SD)",
    color = "Potential Type"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Salvataggio del plot
ggsave(
  filename = file.path(output_dir_saturazione, "Curva_Saturazione_ROC_AUC.png"),
  plot = p_saturazione,
  width = 9, 
  height = 5, 
  dpi = 300
)

cat("\n=== Finito! Dati e plot di saturazione salvati in:", output_dir_saturazione, "===\n")

# --- 8. PLOT SEPARATI DELLE 4 CURVE CON ASSE Y BLOCCATO ---

# Funzione aggiornata con asse Y fisso tra 0.5 e 1.0
salva_plot_singolo_fissato <- function(df, nome_dataset, tipo_potenziale, colore_linea) {
  
  # 1. Filtriamo i dati
  df_filtrato <- df %>%
    filter(Dataset == nome_dataset & potential_type == tipo_potenziale)
  
  # 2. Creazione del plot
  p_singolo <- ggplot(df_filtrato, aes(x = num_pdbs, y = mean_auc)) +
    geom_line(color = colore_linea, linewidth = 1) +
    geom_point(size = 2.5, color = colore_linea) +
    geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc), 
                  width = max(df_filtrato$num_pdbs, na.rm = TRUE) * 0.02, 
                  alpha = 0.5, color = "black") + 
    # === LA MODIFICA È QUI ===
    # Fissiamo i limiti visivi dell'asse Y tra 0.5 e 1.0 senza eliminare dati
    coord_cartesian(ylim = c(0.5, 1.0)) +
    # =========================
  theme_minimal() +
    labs(
      title = paste("Saturazione ROC AUC -", nome_dataset, "( Sum", toupper(tipo_potenziale), ")"),
      x = "Numero stimato di PDB (2187 / n_split)",
      y = "Mean ROC AUC (\u00B1 SD)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      panel.grid.minor = element_blank()
    )
  
  # 3. Nome file (aggiunto un prefisso per distinguerli da quelli non fissati)
  nome_file <- paste0("Curva_Saturazione_YFisso_", gsub(" ", "_", nome_dataset), "_", tipo_potenziale, ".png")
  
  # 4. Salvataggio
  ggsave(
    filename = file.path(output_dir_saturazione, nome_file),
    plot = p_singolo,
    width = 6, 
    height = 4, 
    dpi = 300
  )
}

# Generiamo e salviamo le 4 curve con la nuova funzione
salva_plot_singolo_fissato(df_plot, "Whole Interface", "sym", "dodgerblue2")
salva_plot_singolo_fissato(df_plot, "Whole Interface", "asym", "mediumvioletred")
salva_plot_singolo_fissato(df_plot, "Stratified", "sym", "turquoise4") 
salva_plot_singolo_fissato(df_plot, "Stratified", "asym", "deeppink4")

cat("\n=== Finito! I 4 grafici in scala comparabile (Y da 0.5 a 1.0) sono stati salvati in:", output_dir_saturazione, "===\n")


### SE VOGLIAMO AGGIUNGERE IL PUNTO ARTIFICIALE DOVUTO ALLA ROC AUC CON L'INTERO DATASET

# --- 8. PLOT SEPARATI DELLE 4 CURVE CON ASSE Y BLOCCATO E PUNTO DATASET INTERO ---

# Funzione aggiornata con il parametro auc_intero
salva_plot_singolo_fissato <- function(df, nome_dataset, tipo_potenziale, colore_linea, auc_intero) {
  
  # 1. Filtriamo i dati come prima
  df_filtrato <- df %>%
    filter(Dataset == nome_dataset & potential_type == tipo_potenziale)
  
  # 2. Creiamo un nuovo dataframe contenente SOLO il punto per il dataset intero
  punto_intero <- data.frame(
    mean_auc = auc_intero,
    sd_auc = NA,         # Nessuna deviazione standard per il dataset intero
    num_pdbs = 2187,     # Hardcoded: tutti i PDB
    Dataset = nome_dataset,
    potential_type = tipo_potenziale
  )
  
  # Uniamo il punto finale ai punti degli split e li ordiniamo per asse X
  df_plot_final <- bind_rows(df_filtrato, punto_intero) %>%
    arrange(num_pdbs)
  
  # 3. Creazione del plot
  p_singolo <- ggplot(df_plot_final, aes(x = num_pdbs, y = mean_auc)) +
    # Disegna la linea che ora arriverà fino a x = 2187
    geom_line(color = colore_linea, linewidth = 1) +
    # Punti standard degli split
    geom_point(size = 2.5, color = colore_linea) +
    # Punto speciale per il dataset intero (rombo nero per distinguerlo)
    geom_point(data = punto_intero, aes(x = num_pdbs, y = mean_auc), 
               size = 4, shape = 18, color = "black") +
    # Barre di errore (filtrate per evitare warning dove sd_auc è NA)
    geom_errorbar(data = df_plot_final %>% filter(!is.na(sd_auc)),
                  aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc), 
                  width = 2187 * 0.02, 
                  alpha = 0.5, color = "black") + 
    coord_cartesian(ylim = c(0.5, 1.0)) +
    theme_minimal() +
    labs(
      title = paste("Saturazione ROC AUC -", nome_dataset, "( Sum", toupper(tipo_potenziale), ")"),
      x = "Numero stimato di PDB (2187 / n_split)",
      y = "Mean ROC AUC (\u00B1 SD)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      panel.grid.minor = element_blank()
    )
  
  # 4. Nome file 
  nome_file <- paste0("Curva_Saturazione_YFisso_", gsub(" ", "_", nome_dataset), "_", tipo_potenziale, ".png")
  
  # 5. Salvataggio
  ggsave(
    filename = file.path(output_dir_saturazione, nome_file),
    plot = p_singolo,
    width = 6, 
    height = 4, 
    dpi = 300
  )
}

# Generiamo e salviamo le 4 curve, passando i valori hardcoded dell'intero dataset
# !!! INSERISCI I TUOI VERI VALORI DI AUC AL POSTO DI QUELLI DI ESEMPIO (es. 0.75, 0.82...) !!!
salva_plot_singolo_fissato(df_plot, "Whole Interface", "sym", "dodgerblue2", auc_intero = 0.55)
salva_plot_singolo_fissato(df_plot, "Whole Interface", "asym", "mediumvioletred", auc_intero = 0.70)
salva_plot_singolo_fissato(df_plot, "Stratified", "sym", "turquoise4", auc_intero = 0.60) 
salva_plot_singolo_fissato(df_plot, "Stratified", "asym", "deeppink4", auc_intero = 0.74)

cat("\n=== Finito! I 4 grafici con il punto finale (2187 PDB) sono stati salvati in:", output_dir_saturazione, "===\n")


# --- 9. CONTROLLO OVERFITTING/BIAS: DISTRIBUZIONE DEI PUNTEGGI ---

# Poiché merged_df era intrappolato nella funzione, carichiamo i dati al volo 
# per il dataset "Whole Interface" come campione per il controllo bias.
dockq_raw <- read.csv(path_dockq_af3)
pot_raw_whole <- read.csv(path_pot_af3_whole)
#pot_raw_strat <- read.csv(path_pot_af3_strat)

merged_check <- dockq_raw %>%
  inner_join(pot_raw_whole, by = c("Model" = "pdb"))

# Prepariamo i dati estraendo le due classi estreme
df_density <- merged_check %>%
  filter(DockQ <= 0.24 | DockQ >= 0.81) %>%
  mutate(
    Classe = ifelse(DockQ >= 0.81, "Native-like (High Quality)", "Decoy (Incorrect)")
  )

# Creiamo il grafico a densità
p_density <- ggplot(df_density, aes(x = sum_potential, fill = Classe)) +
  geom_density(alpha = 0.6, color = "black", linewidth = 0.5) +
  facet_wrap(~ potential_type, scales = "free") +
  theme_minimal() +
  scale_fill_manual(values = c("Native-like (High Quality)" = "forestgreen", 
                               "Decoy (Incorrect)" = "firebrick")) +
  labs(
    title = "Controllo Bias: Distribuzione dei Punteggi (Whole-interface)",
    subtitle = "Se le due distribuzioni sono completamente separate, il task è troppo facile.",
    x = "Sum Potential (Score Grezzo)",
    y = "Densità",
    fill = "Tipo di Posa"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Salviamo il plot di controllo
ggsave(
  filename = file.path(output_dir_saturazione, "Controllo_Distribuzione_Punteggi_Whole.png"),
  plot = p_density,
  width = 8, 
  height = 5, 
  dpi = 300
)

cat("\n=== Plot di controllo della distribuzione salvato ===\n")

# Creiamo l'istogramma con i conteggi assoluti
p_histogram <- ggplot(df_density, aes(x = sum_potential, fill = Classe)) +
  # Usiamo bins = 60 (puoi aggiustarlo) e position = "identity" per sovrapporli
  geom_histogram(position = "identity", alpha = 0.6, color = "black", bins = 60, linewidth = 0.2) +
  facet_wrap(~ potential_type, scales = "free") +
  theme_minimal() +
  scale_fill_manual(values = c("Native-like (High Quality)" = "forestgreen", 
                               "Decoy (Incorrect)" = "firebrick")) +
  labs(
    title = "Controllo Bias: Conteggi Assoluti dei Punteggi (Istogramma)",
    subtitle = "L'asse Y mostra il numero reale di pose, rivelando lo sbilanciamento delle classi.",
    x = "Sum Potential (Score Grezzo)",
    y = "Numero di Pose (Count)",
    fill = "Tipo di Posa"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Salviamo il plot
ggsave(
  filename = file.path(output_dir_saturazione, "Controllo_Istogramma_Punteggi.png"),
  plot = p_histogram,
  width = 8, 
  height = 5, 
  dpi = 300
)

cat("\n=== Plot a istogramma salvato! ===\n")

# --- 10. SMASCHERIAMO IL PARADOSSO: ISTOGRAMMA DEGLI Z-SCORE ---

# 1. Calcoliamo lo Z-score usando TUTTE le pose (prima di filtrare per DockQ)
# Raggruppiamo per 'group_num' (il singolo complesso proteico) e 'potential_type'
df_zscore <- merged_check %>%
  group_by(group_num, potential_type) %>%
  mutate(
    # Lo Z-score: (Valore - Media) / Deviazione Standard
    # Ci dice di quante deviazioni standard una posa è migliore o peggiore 
    # rispetto alla media delle pose per QUELLA specifica proteina.
    z_score_potential = (sum_potential - mean(sum_potential, na.rm = TRUE)) / sd(sum_potential, na.rm = TRUE)
  ) %>%
  ungroup()

# 2. Ora filtriamo solo le classi estreme per la visualizzazione
df_zplot <- df_zscore %>%
  filter(DockQ <= 0.24 | DockQ >= 0.81) %>%
  mutate(
    Classe = ifelse(DockQ >= 0.81, "Native-like (High Quality)", "Decoy (Incorrect)")
  )

# 3. Creiamo l'istogramma degli Z-score
p_zscore <- ggplot(df_zplot, aes(x = z_score_potential, fill = Classe)) +
  geom_histogram(position = "identity", alpha = 0.6, color = "black", bins = 60, linewidth = 0.2) +
  facet_wrap(~ potential_type, scales = "free_y") + 
  theme_minimal() +
  scale_fill_manual(values = c("Native-like (High Quality)" = "forestgreen", 
                               "Decoy (Incorrect)" = "firebrick")) +
  labs(
    title = "Verifica del Ranking Locale: Istogramma degli Z-Score",
    subtitle = "Punteggi normalizzati per singola proteina. Annulla il bias della dimensione.",
    x = "Z-Score (0 = Media della proteina, <0 = Migliore della media)",
    y = "Numero di Pose (Count)",
    fill = "Tipo di Posa"
  ) +
  # Fissiamo l'asse X tra -5 e +5 per escludere rarissimi outlier che rovinerebbero la visualizzazione
  coord_cartesian(xlim = c(-5, 5)) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# 4. Salvataggio
ggsave(
  filename = file.path(output_dir_saturazione, "Controllo_Istogramma_ZScore_whole.png"),
  plot = p_zscore,
  width = 8, 
  height = 5, 
  dpi = 300
)

cat("\n=== Plot a istogramma con Z-Score salvato! ===\n")