library(stats)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggdark)

######## PARTE EXTRA PER PLOT SCURO 

dark_minimal <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      # Background
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      
      # Testi
      text = element_text(color = "white"),
      plot.title = element_text(color = "white", hjust = 0.5),
      
      # Assi
      axis.text = element_text(color = "white"),
      axis.title = element_text(color = "white"),
      
      # Griglia
      panel.grid.major = element_line(color = "grey30"),
      panel.grid.minor = element_line(color = "grey20"),
      
      # Legenda
      legend.background = element_rect(fill = "black"),
      legend.key = element_rect(fill = "black"),
      legend.text = element_text(color = "white")
    )
}

##############

path_dockq_hdock = "/Users/lorenzosisti/Downloads/DockQ_HDOCK/DockQ_results_HDOCK.csv"
path_dockq_af3 = "/Users/lorenzosisti/Downloads/DockQ_results_AF3_12_06.csv"

hdock_dockq <- read.csv(path_dockq_hdock)
af3_dockq <- read.csv(path_dockq_af3)
density_hdock <- density(hdock_dockq$DockQ)

hdock_decoy <- nrow(hdock_dockq[hdock_dockq$DockQ <= 0.24,])
hdock_near_native <- nrow(hdock_dockq[hdock_dockq$DockQ >= 0.81,])
hdock_intermediate <- nrow(hdock_dockq) - hdock_decoy - hdock_near_native

# 1. Calcolo del totale
total_poses <- nrow(hdock_dockq)

# 2. Calcolo delle frazioni
frac_decoy <- hdock_decoy / total_poses
frac_near_native <- hdock_near_native / total_poses
frac_intermediate <- hdock_intermediate / total_poses

# 3. Creazione del dataframe per il plot
df_summary <- data.frame(
  Categoria = factor(c("Decoy", "Intermediate", "Near-native"), 
                     levels = c("Decoy", "Intermediate", "Near-native")),
  Frazione = c(frac_decoy, frac_intermediate, frac_near_native),
  Conteggio = c(hdock_decoy, hdock_intermediate, hdock_near_native) # <--- Nuova colonna
)

comp_hdock <- ggplot(df_summary, aes(x = Categoria, y = Frazione)) +
  # Le barre
  geom_col(fill = "dodgerblue2", color = "black", width = 0.6) + 
  
  # I numeri sopra le barre
  geom_text(aes(label = Conteggio), vjust = -0.5, size = 5) + 
  
  # Impostiamo i limiti
  scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
  
  labs(
    y = "Composition fraction",
    x = NULL  # <--- Rimuove il titolo dell'asse X ("Categoria")
  ) +
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    # <--- Modifica qui sotto per le scritte sotto le barre
    axis.text.x = element_text(size = 12, color = "black"),
    # Opzionale: ingrandisce anche il titolo dell'asse Y per coerenza
    axis.title.y = element_text(size = 13) 
  )

# 2. Plotta la PDF stimata (usando geom_density)
pdf_hdock <- ggplot(hdock_dockq, aes(x = DockQ)) +
  geom_density(fill = "dodgerblue2", alpha = 0.6) + 
  labs(
    title = "DockQ PDF - HDOCK",
    x = "DockQ",
    y = "Probability Density Function (PDF)"
  ) +
  theme_minimal() +
  theme(
    # --- TITOLO CENTRATO ---
    plot.title = element_text(hjust = 0.5, size = 12)
  )

p_final <- pdf_hdock + inset_element(
  comp_hdock, 
  left = 0.25,   # Inizia al 60% dell'asse X (verso destra)
  bottom = 0.4, # Inizia al 60% dell'asse Y (verso l'alto)
  right = 1,    # Finisce al margine destro
  top = 1       # Finisce al margine superiore
)

pdf_hdock_dark <- ggplot(hdock_dockq, aes(x = DockQ)) +
  geom_density(fill = "dodgerblue2", alpha = 0.6) + 
  labs(
    title = "DockQ PDF - HDOCK",
    x = "DockQ",
    y = "Probability Density Function (PDF)"
  ) +
  dark_minimal()

comp_hdock_dark <- comp_hdock + dark_minimal(base_size = 9)

p_hdock_dark <- pdf_hdock_dark +
  inset_element(
    comp_hdock_dark,
    left = 0.26,
    bottom = 0.41,
    right = 1,
    top = 1
  )

af3_decoy <- nrow(af3_dockq[af3_dockq$DockQ <= 0.24,])
af3_near_native <- nrow(af3_dockq[af3_dockq$DockQ >= 0.81,])
af3_intermediate <- nrow(af3_dockq) - af3_decoy - af3_near_native

# 1. Calcolo del totale
total_poses <- nrow(af3_dockq)

# 2. Calcolo delle frazioni
frac_decoy <- af3_decoy / total_poses
frac_near_native <- af3_near_native / total_poses
frac_intermediate <- af3_intermediate / total_poses

# 3. Creazione del dataframe per il plot
df_summary <- data.frame(
  Categoria = factor(c("Decoy", "Intermediate", "Near-native"), 
                     levels = c("Decoy", "Intermediate", "Near-native")),
  Frazione = c(frac_decoy, frac_intermediate, frac_near_native),
  Conteggio = c(af3_decoy, af3_intermediate, af3_near_native) # <--- Nuova colonna
)

comp_af3 <- ggplot(df_summary, aes(x = Categoria, y = Frazione)) +
  # Le barre
  geom_col(fill = "dodgerblue2", color = "black", width = 0.6) + 
  
  # I numeri sopra le barre
  geom_text(aes(label = Conteggio), vjust = -0.5, size = 5) + 
  
  # Impostiamo i limiti
  scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
  
  labs(
    y = "Composition fraction",
    x = NULL  # <--- Rimuove il titolo dell'asse X ("Categoria")
  ) +
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    # <--- Modifica qui sotto per le scritte sotto le barre
    axis.text.x = element_text(size = 12, color = "black"),
    # Opzionale: ingrandisce anche il titolo dell'asse Y per coerenza
    axis.title.y = element_text(size = 13) 
  )


# 2. Plotta la PDF stimata (usando geom_density)
pdf_af3 <- ggplot(af3_dockq, aes(x = DockQ)) +
  geom_density(fill = "dodgerblue2", alpha = 0.6) + 
  labs(
    title = "DockQ PDF - AlphaFold 3",
    x = "DockQ",
    y = "Probability Density Function (PDF)"
  ) +
  theme_minimal() +
  theme(
    # --- TITOLO CENTRATO ---
    plot.title = element_text(hjust = 0.5, size = 12)
  )

p_final <- pdf_af3 + inset_element(
  comp_af3, 
  left = 0.25,   # Inizia al 60% dell'asse X (verso destra)
  bottom = 0.4, # Inizia al 60% dell'asse Y (verso l'alto)
  right = 1,    # Finisce al margine destro
  top = 1       # Finisce al margine superiore
)

pdf_af3_dark <- ggplot(af3_dockq, aes(x = DockQ)) +
  geom_density(fill = "dodgerblue2", alpha = 0.6) + 
  labs(
    title = "DockQ PDF - AlphaFold 3",
    x = "DockQ",
    y = "Probability Density Function (PDF)"
  ) +
  dark_minimal()

comp_af3_dark <- comp_af3 + dark_minimal(base_size = 9)

p_af3_dark <- pdf_af3_dark +
  inset_element(
    comp_af3_dark,
    left = 0.26,
    bottom = 0.41,
    right = 1,
    top = 1
  )

