# =========================================================================
#  AUC vs soglia di fnat  —  AF3 e HDOCK, potenziali whole / gr / cdr
# =========================================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(pROC)

# --- 1. IMPOSTAZIONI GLOBALI ---------------------------------------------

# Define the path to a custom function files
source("/Users/lorenzosisti/Documents/Script_ottimizzati_funzioni/functions.R")

output_dir <- "/Users/lorenzosisti/Downloads/AUC_vs_fnat_17_09_024"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

plot_width  <- 9
plot_height <- 7
plot_dpi    <- 300
theme_custom <- theme_minimal() + theme(plot.title = element_blank())

soglia_dockq   <- 0.24                  # DockQ <= soglia  ->  classe positiva (decoy)
soglie_fnat <- seq(0, 1, by = 0.02)
min_per_classe <- 1                    # n minimo di decoy e non-decoy per calcolare l'AUC

# --- 2. DEFINIZIONE DEI DATASET ------------------------------------------
# Per ogni metodo di docking: file DockQ + i tre file di potenziali.
# 'join_col' serve perche' il file CDR usa 'pdb_filename' invece di 'pdb'.

datasets <- list(
  AF3 = list(
    dockq = "/Users/lorenzosisti/Downloads/DockQ_results_AF3_12_06.csv",
    potenziali = list(
      Whole      = list(path = "/Users/lorenzosisti/Downloads/potenziali_statistici_30_06_50_pose/potenziali_whole_per_posa.csv",  join_col = "pdb"),
      Stratified = list(path = "/Users/lorenzosisti/Downloads/potenziali_statistici_30_06_50_pose/potenziali_gr_per_posa.csv",     join_col = "pdb"),
      CDR        = list(path = "/Users/lorenzosisti/Downloads/potenziali_statistici_30_06_50_pose/punteggi_cdr_per_posa.csv",      join_col = "pdb_filename")
    )
  ),
  HDOCK = list(
    dockq = "/Users/lorenzosisti/Downloads/DockQ_HDOCK/DockQ_results_HDOCK.csv",
    potenziali = list(
      Whole      = list(path = "/Users/lorenzosisti/Downloads/potenziali_statistici_30_06_hdock/potenziali_whole_per_posa.csv", join_col = "pdb"),
      Stratified = list(path = "/Users/lorenzosisti/Downloads/potenziali_statistici_30_06_hdock/potenziali_gr_per_posa.csv",    join_col = "pdb"),
      CDR        = list(path = "/Users/lorenzosisti/Downloads/potenziali_statistici_30_06_hdock/punteggi_cdr_per_posa.csv",     join_col = "pdb_filename")
    )
  )
)

# --- 3. FUNZIONI ----------------------------------------------------------

# 3a. Unisce DockQ e potenziali e aggiunge la classe vera.
prepara_dati <- function(path_dockq, path_potenziale, join_col) {
  
  dockq_scores     <- read.csv(path_dockq)
  potential_scores <- read.csv(path_potenziale)
  
  by_vec <- setNames(join_col, "Model")   # es. c("Model" = "pdb_filename")
  
  dockq_scores %>%
    left_join(potential_scores, by = by_vec) %>%
    mutate(true_class = ifelse(DockQ <= soglia_dockq, 1, 0))
}


# 3b. Stima la direzione della ROC UNA VOLTA sul dataset completo,
#     cosi' tutte le soglie usano la stessa convenzione e sono confrontabili.
stima_direzione <- function(df, score_col) {
  r <- roc(df$true_class, df[[score_col]], quiet = TRUE)
  r$direction
}

# 3c. Cuore del calcolo: per ogni soglia di fnat, sottoinsieme + AUC.
#     Restituisce NA quando una delle due classi e' troppo poco popolata.
auc_vs_fnat <- function(df, score_col, direzione, soglie = soglie_fnat) {
  
  df <- df[!is.na(df$fnat) & !is.na(df[[score_col]]) & !is.na(df$true_class), ]
  
  righe <- lapply(soglie, function(soglia) {
    
    sub   <- df[df$fnat > soglia, ]
    n_pos <- sum(sub$true_class == 1)   # decoy
    n_neg <- sum(sub$true_class == 0)   # non-decoy
    
    # caso degenere: non ha senso calcolare una ROC
    if (n_pos < min_per_classe || n_neg < min_per_classe) {
      return(data.frame(soglia_fnat = soglia, n_tot = nrow(sub),
                        n_pos = n_pos, n_neg = n_neg,
                        auc = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
    }
    
    r  <- roc(sub$true_class, sub[[score_col]],
              direction = direzione, quiet = TRUE)
    ci <- as.numeric(ci.auc(r, method = "delong"))   # [inf, stima, sup]
    
    data.frame(soglia_fnat = soglia, n_tot = nrow(sub),
               n_pos = n_pos, n_neg = n_neg,
               auc = as.numeric(auc(r)), ci_low = ci[1], ci_high = ci[3])
  })
  
  bind_rows(righe)
}

# 3d. Grafico: pannello superiore = AUC, pannello inferiore = numerosita'.
plot_auc_vs_fnat <- function(risultati, etichetta) {
  
  p_auc <- ggplot(risultati, aes(x = soglia_fnat, y = auc, color = potenziale)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = potenziale),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = NULL, y = "ROC AUC", color = "Potenziale", fill = "Potenziale") +
    theme_custom +
    theme(legend.position = "right")
  
  p_n <- ggplot(risultati, aes(x = soglia_fnat)) +
    geom_line(aes(y = n_pos, linetype = "Decoy (DockQ <= 0.24)"), linewidth = 0.7) +
    geom_line(aes(y = n_neg, linetype = "Non decoy"),             linewidth = 0.7) +
    labs(x = paste0("Soglia di fnat (pose con fnat > soglia) — ", etichetta),
         y = "N pose", linetype = NULL) +
    theme_custom +
    theme(legend.position = "right")
  
  p_auc / p_n + plot_layout(heights = c(3, 1))
}

# --- 4. CICLO SU METODI E POTENZIALI --------------------------------------

risultati_completi <- list()

for (metodo in names(datasets)) {
  
  cfg <- datasets[[metodo]]
  
  # se il file DockQ non e' disponibile, salto il metodo senza bloccare lo script
  if (!file.exists(cfg$dockq)) {
    message("File DockQ non trovato per ", metodo, ": salto.")
    next
  }
  
  for (nome_pot in names(cfg$potenziali)) {
    
    pot <- cfg$potenziali[[nome_pot]]
    if (!file.exists(pot$path)) {
      message("File potenziale non trovato: ", pot$path, " — salto.")
      next
    }
    
    df <- prepara_dati(cfg$dockq, pot$path, pot$join_col)
    
    # due score per ogni potenziale: sym e asym
    res_pot <- bind_rows(lapply(c("mean_sym", "mean_asym"), function(sc) {
      dir_sc <- stima_direzione(df, sc)
      auc_vs_fnat(df, sc, dir_sc) %>%
        mutate(metodo     = metodo,
               potenziale = paste0(nome_pot, " - ",
                                   ifelse(sc == "mean_sym", "Sym", "Asym")),
               score      = sc,
               direzione  = dir_sc)
    }))
    
    risultati_completi[[paste(metodo, nome_pot, sep = "_")]] <- res_pot
    
    # --- grafico singolo (uno per metodo x potenziale) ---
    p <- plot_auc_vs_fnat(res_pot, paste(metodo, nome_pot))
    print(p)
    
    ggsave(
      filename = file.path(output_dir,
                           paste0(metodo, "_AUCvsFnat_", nome_pot, ".png")),
      plot = p, width = plot_width, height = plot_height, dpi = plot_dpi
    )
  }
}

# tabella unica con tutti i punti, utile per controlli e per l'export
risultati_df <- bind_rows(risultati_completi)
write.csv(risultati_df,
          file.path(output_dir, "AUC_vs_fnat_tutti.csv"), row.names = FALSE)

# --- 5. GRAFICO RIASSUNTIVO (facoltativo) ---------------------------------
# Tutte le combinazioni affiancate, un pannello per metodo di docking.

p_riassunto <- ggplot(risultati_df,
                      aes(x = soglia_fnat, y = auc, color = potenziale)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  facet_wrap(~ metodo) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Soglia di fnat (pose con fnat > soglia)",
       y = "ROC AUC", color = "Potenziale") +
  theme_custom +
  theme(legend.position = "right")

print(p_riassunto)

ggsave(file.path(output_dir, "AUCvsFnat_riassunto.png"),
       p_riassunto, width = plot_width + 3, height = plot_height,
       dpi = plot_dpi)

# =========================================================================
#  APPENDICE — AUC vs soglia su due metriche: fnat e %fw
#  Da eseguire DOPO lo script principale (riusa theme_custom, soglia_dockq,
#  min_per_classe, stima_direzione, output_dir).
# =========================================================================

library(data.table)

# --- A1. PARAMETRI SPECIFICI DELL'APPENDICE ------------------------------

min_contatti <- 1   # n minimo di contatti per considerare affidabile la %fw

# Griglie di soglie. ATTENZIONE alla lettura dell'asse x:
#   - fnat:   ci si sposta a destra = filtro PIU' STRINGENTE (tengo fnat > soglia)
#   - pct_fw: ci si sposta a destra = filtro PIU' PERMISSIVO (tengo %fw < soglia)
metriche <- list(
  fnat = list(
    col    = "fnat",
    verso  = "greater",
    soglie = seq(0, 0.60, by = 0.02),
    label  = "Soglia di fnat  (tengo le pose con fnat > soglia)"
  ),
  pct_fw = list(
    col    = "pct_fw",
    verso  = "less",
    soglie = seq(5, 100, by = 5),
    label  = "Soglia di %fw  (tengo le pose con %fw < soglia)"
  )
)

# File dei contatti prodotti dallo script di filtraggio CDR, uno per metodo.
# NB: per HDOCK va rieseguito quello script puntando alla directory HDOCK.
path_contacts <- list(
  AF3   = "/Users/lorenzosisti/Downloads/cdrs_filtered_docked_structure/df_contacts.csv",
  HDOCK = "/Users/lorenzosisti/Downloads/hdock_cdrs_filtered_docked_structure/df_contacts.csv"
)

# --- A2. DA df_contacts A UNA TABELLA %fw PER POSA ------------------------

calcola_pct_fw <- function(path_contacts_csv) {
  
  dt <- fread(path_contacts_csv)
  
  # Due definizioni affiancate:
  #  - pct_fw     : % di COPPIE di contatto che coinvolgono un residuo fw
  #  - pct_fw_res : % di RESIDUI Ab in contatto che sono fw (meno sensibile
  #                 al fatto che un singolo residuo tocchi molti residui Ag)
  res <- dt[, .(
    n_total     = .N,
    n_fw        = sum(region_ab == "fw"),
    n_res_ab    = uniqueN(paste(chain_ab, resno_ab)),
    n_res_ab_fw = uniqueN(paste(chain_ab, resno_ab)[region_ab == "fw"])
  ), by = pdb_id][
    , `:=`(pct_fw     = 100 * n_fw / n_total,
           pct_fw_res = 100 * n_res_ab_fw / n_res_ab)
  ]
  
  # pdb_id = "<nome_file>_<H>_<L>_<Ag>"  ->  Model = "<nome_file>.pdb"
  res[, Model := paste0(sub("_[^_]+_[^_]+_[^_]+$", "", pdb_id), ".pdb")]
  
  res[, .(Model, pdb_id, n_total, pct_fw, pct_fw_res)]
}

# --- A3. PREPARAZIONE DATI (versione estesa, con %fw) --------------------

prepara_dati_esteso <- function(path_dockq, path_potenziale, join_col,
                                path_contacts_csv) {
  
  dockq_scores     <- read.csv(path_dockq)
  potential_scores <- read.csv(path_potenziale)
  fw_scores        <- calcola_pct_fw(path_contacts_csv)
  
  by_vec <- setNames(join_col, "Model")
  
  df <- dockq_scores %>%
    left_join(potential_scores, by = by_vec) %>%
    left_join(as.data.frame(fw_scores), by = "Model") %>%
    mutate(true_class = ifelse(DockQ <= soglia_dockq, 1, 0))
  
  # diagnostica: quante pose non hanno trovato la controparte nei contatti
  n_mancanti <- sum(is.na(df$pct_fw))
  if (n_mancanti > 0) {
    message("  Pose senza %fw (nessun match in df_contacts): ", n_mancanti)
  }
  
  # pose con interfaccia troppo piccola: %fw non affidabile -> la annullo
  df$pct_fw[!is.na(df$n_total) & df$n_total < min_contatti] <- NA_real_
  
  df
}

# --- A4. CALCOLO AUC AL VARIARE DI UNA METRICA GENERICA ------------------

auc_vs_metrica <- function(df, score_col, direzione,
                           metrica_col, verso, soglie) {
  
  df <- df[!is.na(df[[metrica_col]]) &
             !is.na(df[[score_col]])   &
             !is.na(df$true_class), ]
  
  righe <- lapply(soglie, function(soglia) {
    
    # unica differenza tra le due metriche: il verso della disuguaglianza
    tieni <- if (verso == "greater") df[[metrica_col]] >  soglia
    else                    df[[metrica_col]] <  soglia
    sub   <- df[tieni, ]
    
    n_pos <- sum(sub$true_class == 1)   # decoy
    n_neg <- sum(sub$true_class == 0)   # non decoy
    
    if (n_pos < min_per_classe || n_neg < min_per_classe) {
      return(data.frame(soglia = soglia, n_tot = nrow(sub),
                        n_pos = n_pos, n_neg = n_neg,
                        auc = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
    }
    
    r  <- roc(sub$true_class, sub[[score_col]],
              direction = direzione, quiet = TRUE)
    ci <- as.numeric(ci.auc(r, method = "delong"))
    
    data.frame(soglia = soglia, n_tot = nrow(sub),
               n_pos = n_pos, n_neg = n_neg,
               auc = as.numeric(auc(r)), ci_low = ci[1], ci_high = ci[3])
  })
  
  bind_rows(righe)
}

# --- A5. GRAFICO (stessa struttura di prima, asse x parametrico) ---------

plot_auc_vs_metrica <- function(risultati, x_label, etichetta) {
  
  p_auc <- ggplot(risultati, aes(x = soglia, y = auc, color = potenziale)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
    geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = potenziale),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = NULL, y = "ROC AUC", color = "Potenziale", fill = "Potenziale") +
    theme_custom +
    theme(legend.position = "right")
  
  p_n <- ggplot(risultati, aes(x = soglia)) +
    geom_line(aes(y = n_pos, linetype = "Decoy (DockQ <= 0.24)"), linewidth = 0.7) +
    geom_line(aes(y = n_neg, linetype = "Non decoy"),             linewidth = 0.7) +
    labs(x = paste0(x_label, "  —  ", etichetta), y = "N pose", linetype = NULL) +
    theme_custom +
    theme(legend.position = "right")
  
  p_auc / p_n + plot_layout(heights = c(3, 1))
}

# --- A6. CICLO SU METODI x POTENZIALI x METRICHE -------------------------

risultati_appendice <- list()

for (metodo in names(datasets)) {
  
  cfg <- datasets[[metodo]]
  
  if (!file.exists(cfg$dockq) || !file.exists(path_contacts[[metodo]])) {
    message("Dati incompleti per ", metodo, " (DockQ o df_contacts): salto.")
    next
  }
  
  for (nome_pot in names(cfg$potenziali)) {
    
    pot <- cfg$potenziali[[nome_pot]]
    if (!file.exists(pot$path)) {
      message("File potenziale non trovato: ", pot$path, " — salto.")
      next
    }
    
    message("Elaboro: ", metodo, " / ", nome_pot)
    df <- prepara_dati_esteso(cfg$dockq, pot$path, pot$join_col,
                              path_contacts[[metodo]])
    
    for (nome_metrica in names(metriche)) {
      
      m <- metriche[[nome_metrica]]
      
      res <- bind_rows(lapply(c("mean_sym", "mean_asym"), function(sc) {
        dir_sc <- stima_direzione(df, sc)     # direzione fissata sul dataset pieno
        auc_vs_metrica(df, sc, dir_sc, m$col, m$verso, m$soglie) %>%
          mutate(metodo     = metodo,
                 metrica    = nome_metrica,
                 potenziale = paste0(nome_pot, " - ",
                                     ifelse(sc == "mean_sym", "Sym", "Asym")))
      }))
      
      risultati_appendice[[paste(metodo, nome_pot, nome_metrica, sep = "_")]] <- res
      
      p <- plot_auc_vs_metrica(res, m$label, paste(metodo, nome_pot))
      print(p)
      
      ggsave(
        filename = file.path(output_dir,
                             paste0(metodo, "_AUCvs", nome_metrica,
                                    "_", nome_pot, ".png")),
        plot = p, width = plot_width, height = plot_height, dpi = plot_dpi
      )
    }
  }
}

risultati_appendice_df <- bind_rows(risultati_appendice)
write.csv(risultati_appendice_df,
          file.path(output_dir, "AUC_vs_metriche_tutti.csv"), row.names = FALSE)

file.exists(datasets$HDOCK$dockq)
file.exists(path_contacts$HDOCK)

# --- A7. RIASSUNTO: metodo x metrica --------------------------------------

p_riassunto2 <- ggplot(risultati_appendice_df,
                       aes(x = soglia, y = auc, color = potenziale)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.4) +
  facet_grid(metodo ~ metrica, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Soglia sulla metrica di filtro", y = "ROC AUC",
       color = "Potenziale") +
  theme_custom +
  theme(legend.position = "right")

print(p_riassunto2)

ggsave(file.path(output_dir, "AUCvsMetriche_riassunto.png"),
       p_riassunto2, width = plot_width + 4, height = plot_height + 2,
       dpi = plot_dpi)

summary(dockq_scores$fnat)
quantile(dockq_scores$fnat, probs = seq(0, 1, 0.1), na.rm = TRUE)

# Ci sono decoy con fnat alto? Se sì, quali e quanti?
decoy_fnat_alto <- dockq_scores[dockq_scores$DockQ <= 0.24 & dockq_scores$fnat > 0.25, ]
nrow(decoy_fnat_alto)
head(decoy_fnat_alto)


# =========================================================================
#  APPENDICE B — Test di DeLong appaiati, su soglia di %fw
#  B1: Sym vs Asym (stesso potenziale, stesso metodo)
#  B2: Asym — CDR vs Whole e Stratified vs Whole (stesso metodo)
#  Richiede A2/A3 (calcola_pct_fw, prepara_dati_esteso) gia' definite.
# =========================================================================

# --- B0. Funzione di supporto: DeLong appaiato ad una soglia -------------
delong_a_soglia <- function(df, filtro_col, soglia, score_A, score_B,
                            dir_A, dir_B, verso = "less") {
  
  tieni <- !is.na(df[[filtro_col]]) &
             (if (verso == "greater") df[[filtro_col]] > soglia else df[[filtro_col]] < soglia)
  sub   <- df[tieni & !is.na(df[[score_A]]) & !is.na(df[[score_B]]) & !is.na(df$true_class), ]
  
  n_pos <- sum(sub$true_class == 1, na.rm = TRUE)
  n_neg <- sum(sub$true_class == 0, na.rm = TRUE)
  
  if (n_pos < min_per_classe || n_neg < min_per_classe) {
    return(data.frame(soglia = soglia, n_tot = nrow(sub), n_pos = n_pos, n_neg = n_neg,
                      auc_A = NA_real_, auc_B = NA_real_,
                      diff_auc = NA_real_, p_value = NA_real_))
  }
  
  r_A <- roc(sub$true_class, sub[[score_A]], direction = dir_A, quiet = TRUE)
  r_B <- roc(sub$true_class, sub[[score_B]], direction = dir_B, quiet = TRUE)
  
  # tryCatch: se roc.test fallisce (es. varianza degenere con AUC=1 su entrambe
  # le curve), non blocchiamo il ciclo ma registriamo NA e continuiamo
  test_result <- tryCatch({
    t <- suppressWarnings(roc.test(r_A, r_B, method = "delong", paired = TRUE))
    t$p.value
  }, error = function(e) {
    message("  DeLong non calcolabile a soglia ", soglia,
            " (n_pos=", n_pos, ", n_neg=", n_neg, "): ", conditionMessage(e))
    NA_real_
  })
  
  data.frame(soglia = soglia, n_tot = nrow(sub), n_pos = n_pos, n_neg = n_neg,
             auc_A = as.numeric(auc(r_A)), auc_B = as.numeric(auc(r_B)),
             diff_auc = as.numeric(auc(r_A)) - as.numeric(auc(r_B)),
             p_value = test_result)
}

# griglia di soglie sulla %fw, coerente con quella gia' usata in A1
soglie_pct_fw <- metriche$pct_fw$soglie

# --- B1. SYM vs ASYM, su soglia %fw, per ogni metodo e potenziale -------

risultati_sym_vs_asym <- list()

for (metodo in names(datasets)) {
  
  cfg <- datasets[[metodo]]
  if (!file.exists(cfg$dockq) || !file.exists(path_contacts[[metodo]])) next
  
  for (nome_pot in names(cfg$potenziali)) {
    
    pot <- cfg$potenziali[[nome_pot]]
    if (!file.exists(pot$path)) next
    
    df <- prepara_dati_esteso(cfg$dockq, pot$path, pot$join_col, path_contacts[[metodo]])
    
    dir_sym  <- stima_direzione(df, "mean_sym")
    dir_asym <- stima_direzione(df, "mean_asym")
    
    res <- bind_rows(lapply(soglie_pct_fw, function(s) {
      delong_a_soglia(df, "pct_fw", s, "mean_asym", "mean_sym",
                      dir_asym, dir_sym, verso = "less")
      # A = Asym, B = Sym -> diff_auc > 0 significa Asym migliore di Sym
    })) %>%
      mutate(metodo = metodo, potenziale = nome_pot, confronto = "Asym vs Sym")
    
    risultati_sym_vs_asym[[paste(metodo, nome_pot, sep = "_")]] <- res
  }
}

risultati_sym_vs_asym_df <- bind_rows(risultati_sym_vs_asym) %>%
  group_by(metodo, potenziale) %>%
  mutate(p_value_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup()

write.csv(risultati_sym_vs_asym_df,
          file.path(output_dir, "DeLong_AsymVsSym_pctfw.csv"), row.names = FALSE)

p_sym_asym <- ggplot(risultati_sym_vs_asym_df,
                     aes(x = soglia, y = diff_auc, color = p_value_BH < 0.05)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 1.8) +
  geom_line(aes(group = 1), color = "grey70", linewidth = 0.4) +
  facet_grid(metodo ~ potenziale) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "firebrick"),
                     name = "Significativo\n(BH < 0.05)") +
  labs(x = "Soglia di %fw (pose con %fw < soglia)", y = "AUC(Asym) - AUC(Sym)") +
  theme_custom

print(p_sym_asym)
ggsave(file.path(output_dir, "DeLong_AsymVsSym_pctfw.png"), p_sym_asym,
       width = plot_width + 3, height = plot_height, dpi = plot_dpi)

# --- B2. Solo ASYM: CDR vs Whole e Stratified vs Whole, su soglia %fw ---

risultati_confronto_potenziali <- list()

for (metodo in names(datasets)) {
  
  cfg <- datasets[[metodo]]
  if (!file.exists(cfg$dockq) || !file.exists(path_contacts[[metodo]])) next
  
  pot_whole <- cfg$potenziali[["Whole"]]
  if (is.null(pot_whole) || !file.exists(pot_whole$path)) next
  df_whole <- prepara_dati_esteso(cfg$dockq, pot_whole$path, pot_whole$join_col, path_contacts[[metodo]])
  
  for (nome_pot in c("CDR", "Stratified")) {
    
    pot <- cfg$potenziali[[nome_pot]]
    if (is.null(pot) || !file.exists(pot$path)) next
    
    df_alt <- prepara_dati_esteso(cfg$dockq, pot$path, pot$join_col, path_contacts[[metodo]])
    
    # allinea sulle stesse pose; %fw viene da entrambi i df (stesso df_contacts,
    # quindi identico), lo teniamo da df_whole come riferimento unico
    df_join <- inner_join(
      df_whole %>% select(Model, pct_fw, true_class, mean_asym_whole = mean_asym),
      df_alt   %>% select(Model, mean_asym_alt = mean_asym),
      by = "Model"
    )
    
    dir_whole <- stima_direzione(df_whole, "mean_asym")
    dir_alt   <- stima_direzione(df_alt,   "mean_asym")
    
    res <- bind_rows(lapply(soglie_pct_fw, function(s) {
      delong_a_soglia(df_join, "pct_fw", s, "mean_asym_alt", "mean_asym_whole",
                      dir_alt, dir_whole, verso = "less")
      # A = potenziale alternativo (CDR/Stratified), B = Whole
    })) %>%
      mutate(metodo = metodo, confronto = paste0(nome_pot, " vs Whole (Asym)"))
    
    risultati_confronto_potenziali[[paste(metodo, nome_pot, sep = "_")]] <- res
  }
}

risultati_confronto_potenziali_df <- bind_rows(risultati_confronto_potenziali) %>%
  group_by(metodo, confronto) %>%
  mutate(p_value_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup()

write.csv(risultati_confronto_potenziali_df,
          file.path(output_dir, "DeLong_CDRvsWhole_StratVsWhole_Asym_pctfw.csv"), row.names = FALSE)

p_confronto_pot <- ggplot(risultati_confronto_potenziali_df,
                          aes(x = soglia, y = diff_auc, color = p_value_BH < 0.05)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 1.8) +
  geom_line(aes(group = 1), color = "grey70", linewidth = 0.4) +
  facet_grid(metodo ~ confronto) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "firebrick"),
                     name = "Significativo\n(BH < 0.05)") +
  labs(x = "Soglia di %fw (pose con %fw < soglia)", y = "AUC(alternativa) - AUC(Whole)") +
  theme_custom

print(p_confronto_pot)
ggsave(file.path(output_dir, "DeLong_CDRvsWhole_StratVsWhole_Asym_pctfw.png"), p_confronto_pot,
       width = plot_width + 3, height = plot_height, dpi = plot_dpi)
