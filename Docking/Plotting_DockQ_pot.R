library(dplyr)
library(ggplot2)
library(ggExtra)
library(patchwork)

dockq_scores_hdock <- read.csv("/Users/lorenzosisti/Downloads/DockQ_HDOCK/DockQ_results_HDOCK.csv")
whole_int_pot_scores_hdock <- read.csv("/Users/lorenzosisti/Downloads/potenziali_statistici_HDOCK/potenziali_per_modello.csv")

merged_df <- dockq_scores_hdock %>%
  left_join(whole_int_pot_scores_hdock, by = c("Model" = "pdb"))

p1 <- ggplot(merged_df, aes(x = DockQ, y = sum_sym)) +
  geom_point(alpha = 0.4, size = 1) +   # alpha=trasparenza per gestire 22k punti
  theme_minimal() +
  labs(
    title = "Scatterplot completo: DockQ vs Symmetric Potential",
    x = "DockQ score",
    y = "Symmetric interface potential (sum_sym)"
  )

p2 <- ggplot(merged_df, aes(x = DockQ, y = sum_asym)) +
  geom_point(alpha = 0.4, size = 1) +   # alpha=trasparenza per gestire 22k punti
  theme_minimal() +
  labs(
    title = "Scatterplot completo: DockQ vs Asymmetric Potential",
    x = "DockQ score",
    y = "Symmetric interface potential (sum_asym)"
  )

p1
p2
p1 | p2

p1_marginal <- ggMarginal(p1, type = "density", fill = "gray")
p2_marginal <- ggMarginal(p2, type = "density", fill = "gray")

p1_marginal
p2_marginal

cor(merged_df$DockQ, merged_df$sum_sym)
cor(merged_df$DockQ, merged_df$sum_asym)


###

dockq_scores_AF3 <- read.csv("/Users/lorenzosisti/Downloads/DockQ_HDOCK/DockQ_results_AF3.csv")
whole_int_pot_scores_AF3 <- read.csv("/Users/lorenzosisti/Downloads/potenziali_statistici_AF3/potenziali_per_modello_AF3.csv")

merged_df <- dockq_scores_AF3 %>%
  left_join(whole_int_pot_scores_AF3, by = c("Model" = "pdb"))

p1 <- ggplot(merged_df, aes(x = DockQ, y = sum_sym)) +
  geom_point(alpha = 0.4, size = 1) +   # alpha=trasparenza per gestire 22k punti
  #geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.7) +   # ✅ linea di regressione
  theme_minimal() +
  labs(
    title = "Scatterplot completo: DockQ vs Symmetric Potential",
    x = "DockQ score",
    y = "Symmetric interface potential (sum_sym)"
  )

p2 <- ggplot(merged_df, aes(x = DockQ, y = sum_asym)) +
  geom_point(alpha = 0.4, size = 1) +   # alpha=trasparenza per gestire 22k punti
  #geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.7) +   # ✅ linea di regressione
  theme_minimal() +
  labs(
    title = "Scatterplot completo: DockQ vs Asymmetric Potential",
    x = "DockQ score",
    y = "Symmetric interface potential (sum_asym)"
  )

p1
p2
p1 | p2

p1_marginal <- ggMarginal(p1, type = "density", fill = "gray")
p2_marginal <- ggMarginal(p2, type = "density", fill = "gray")

p1_marginal
p2_marginal

cor(merged_df$DockQ, merged_df$sum_sym)
cor(merged_df$DockQ, merged_df$sum_asym)
