library(data.table)
library(ggplot2)
library(latex2exp)

base    <- "/Users/lorenzosisti/TiNDER/data/repeats"
AUC_COL <- "auc_asym"          # or "auc_sym"
if (!exists("tema_bianco")) tema_bianco <- theme_minimal()

load_rep <- function(root_name, label) {
  rep_dirs <- list.dirs(base, recursive = FALSE, full.names = TRUE)
  rep_dirs <- rep_dirs[grepl("rep_[0-9]+$", rep_dirs)]
  
  rbindlist(lapply(rep_dirs, function(rd) {
    f <- file.path(rd, root_name, paste0("roc_auc_", root_name, ".csv"))
    if (!file.exists(f)) return(NULL)
    d <- fread(f)
    groups_dir <- file.path(rd, "convergence_asymmetric_GR", "groups")
    d[, n_pdb := mapply(function(s, n, g)
      length(readLines(file.path(groups_dir,
                                 sprintf("seed_%s_split_%s_group_%s.txt", s, n, g)))),
      seed, split, group)]
    d[, `:=`(rep = basename(rd), potential = label)]
    d
  }))
}

roc_all <- rbindlist(list(
  load_rep("convergence_asymmetric_GR",  "g(r)"),
  load_rep("convergence_asymmetric_CDR", "CDR")
))

# Step 1: mean over the groups of each repetition
per_rep <- roc_all[, .(rep_mean = mean(get(AUC_COL), na.rm = TRUE),
                       n_pdb    = mean(n_pdb)),
                   by = .(potential, rep, split)]

# # Step 2: mean and sd across repetitions of those means
# aggregati <- per_rep[, .(avg_mean = mean(rep_mean),
#                          avg_sd   = sd(rep_mean),
#                          n_pdb    = mean(n_pdb),
#                          n_rep    = .N),
#                      by = .(potential, split)]
# setorder(aggregati, potential, n_pdb)
# print(aggregati)

# Step 2: mean of the repetition means; sd of the single-group AUCs
sd_single <- roc_all[, .(sd_auc = sd(get(AUC_COL), na.rm = TRUE)),
                     by = .(potential, split)]

aggregati <- per_rep[, .(avg_mean = mean(rep_mean),
                         avg_sd   = sd(rep_mean),     # sd across repetition means (kept)
                         n_pdb    = mean(n_pdb),
                         n_rep    = .N),
                     by = .(potential, split)]
aggregati <- merge(aggregati, sd_single, by = c("potential", "split"))
setorder(aggregati, potential, n_pdb)

make_plot <- function(dat) {
  ggplot(dat, aes(x = n_pdb, y = avg_mean)) +
    geom_line(color = "#CC4678FF", linewidth = 1) +
    geom_point(color = "#CC4678FF", size = 2) +
    geom_errorbar(aes(ymin = avg_mean - avg_sd, ymax = avg_mean + avg_sd),
                  width = diff(range(dat$n_pdb)) * 0.01,
                  color = "#CC4678FF", alpha = 0.6) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
    labs(x = TeX("Number of PDB files used to compute the statistical potentials"),
         y = TeX("Average ROC AUC of the asymmetric potential ($\\pm \\sigma$)")) +
    tema_bianco
}

plot_gr  <- make_plot(aggregati[potential == "g(r)"])
plot_cdr <- make_plot(aggregati[potential == "CDR"])
print(plot_gr); print(plot_cdr)

ggsave(file.path(base, "roc_auc_asym_convergence_GR.png"),  plot_gr,  width = 7, height = 5, dpi = 300, bg = "white")
ggsave(file.path(base, "roc_auc_asym_convergence_CDR.png"), plot_cdr, width = 7, height = 5, dpi = 300, bg = "white")

# ---------------------------------------------------------------
# Scatterplot of the individual ROC AUCs (one point per group x repetition)
# ---------------------------------------------------------------
roc_pts <- copy(roc_all)
roc_pts[, auc := get(AUC_COL)]
roc_pts <- roc_pts[!is.na(auc)]

make_scatter <- function(dat, agg) {
  ggplot() +
    geom_jitter(data = dat, aes(x = n_pdb, y = auc),
                color = "#CC4678FF", alpha = 0.13, size = 1.5,
                width = diff(range(dat$n_pdb)) * 0.005, height = 0) +
    geom_errorbar(data = agg,
                  aes(x = n_pdb, ymin = avg_mean - sd_auc, ymax = avg_mean + sd_auc),
                  color = "#85294CFF", width = 0, linewidth = 1.0, alpha = 0.8) +
    geom_line(data = agg, aes(x = n_pdb, y = avg_mean),
              color = "#85294CFF", linewidth = 0.8) +
    geom_point(data = agg, aes(x = n_pdb, y = avg_mean),
               color = "#85294CFF", size = 2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
    labs(x = "Number of PDB files used to compute the statistical potentials",
         y = "ROC AUC from the regional asymmetric potential") +
    tema_bianco
}

scatter_gr  <- make_scatter(roc_pts[potential == "g(r)"], aggregati[potential == "g(r)"])
scatter_cdr <- make_scatter(roc_pts[potential == "CDR"],  aggregati[potential == "CDR"])
print(scatter_gr); print(scatter_cdr)

ggsave(file.path(base, "roc_auc_asym_scatter_GR.png"),  scatter_gr,  width = 7, height = 5, dpi = 300, bg = "white")
ggsave(file.path(base, "roc_auc_asym_scatter_CDR.png"), scatter_cdr, width = 7, height = 5, dpi = 300, bg = "white")

