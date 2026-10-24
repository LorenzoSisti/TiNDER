library(data.table)
library(pROC)

### ------------------------------------------------------------------
### HEADER: the only part to edit
### ------------------------------------------------------------------
base <- "/Users/lorenzosisti/TiNDER/data"

jobs <- list(
  list(root = file.path(base, "convergence_asymmetric_GR"),  csv_name = "layer_score_dockq.csv"),
  list(root = file.path(base, "convergence_asymmetric_CDR"), csv_name = "cdr_score_dockq.csv")
)

DOCKQ_THR <- 0.81
DIRECTION <- "auto"   # "auto" replicates your original script; "<" is the safer choice (see note)

### ------------------------------------------------------------------
### FUNCTIONS
### ------------------------------------------------------------------
get_auc <- function(df, col) {
  d <- df[!is.na(df[[col]]) & !is.na(df$DockQ), ]
  if (nrow(d) == 0 || length(unique(d$true_class)) < 2) return(NA_real_)
  as.numeric(auc(roc(d$true_class, d[[col]], direction = DIRECTION, quiet = TRUE)))
}

compute_roc_table <- function(root, csv_name) {
  files <- list.files(root, pattern = paste0("^", csv_name, "$"),
                      recursive = TRUE, full.names = TRUE)
  message(basename(root), ": ", length(files), " files found")
  
  rbindlist(lapply(files, function(f) {
    df <- as.data.frame(fread(f))
    df$true_class <- ifelse(df$DockQ <= DOCKQ_THR, 1, 0)
    
    # seed / split / group parsed from the path
    m <- regmatches(f, regexec("seed_([0-9]+)/split_([0-9]+)/group_([0-9]+)", f))[[1]]
    
    data.table(
      seed     = as.integer(m[2]),
      split    = as.integer(m[3]),
      group    = as.integer(m[4]),
      n_poses  = sum(!is.na(df$mean_sym)),
      auc_sym  = get_auc(df, "mean_sym"),
      auc_asym = get_auc(df, "mean_asym"),
      dir      = dirname(f)
    )
  }))
}

### ------------------------------------------------------------------
### RUN
### ------------------------------------------------------------------
roc_tables <- list()

for (j in jobs) {
  res <- compute_roc_table(j$root, j$csv_name)
  setorder(res, split, seed, group)
  
  out_file <- file.path(j$root, paste0("roc_auc_", basename(j$root), ".csv"))
  fwrite(res, out_file)
  message("Saved: ", out_file)
  
  roc_tables[[basename(j$root)]] <- res
  
  print(res[, .(mean_auc_sym  = mean(auc_sym,  na.rm = TRUE),
                mean_auc_asym = mean(auc_asym, na.rm = TRUE),
                n_groups      = .N), by = split])
}

# Tables stay available in the environment:
# roc_tables$convergence_asymmetric_GR
# roc_tables$convergence_asymmetric_CDR

library(data.table)
library(ggplot2)
library(latex2exp)

base       <- "/Users/lorenzosisti/TiNDER/data"
groups_dir <- file.path(base, "convergence_asymmetric_GR", "groups")  # lists made by make_groups.R

if (!exists("tema_bianco")) tema_bianco <- theme_minimal()   # use your own theme if already defined

# ---- Load ROC tables and add the number of PDBs per group ----
load_roc <- function(root_name, label) {
  d <- fread(file.path(base, root_name, paste0("roc_auc_", root_name, ".csv")))
  d[, n_pdb := mapply(function(s, n, g) {
    f <- file.path(groups_dir, sprintf("seed_%s_split_%s_group_%s.txt", s, n, g))
    length(readLines(f))
  }, seed, split, group)]
  d[, potential := label]
  d
}

roc_all <- rbindlist(list(
  load_roc("convergence_asymmetric_GR",  "g(r)"),
  load_roc("convergence_asymmetric_CDR", "CDR")
))

# ---- Aggregate: mean and sd of the asymmetric AUC across seeds/groups ----
aggregati <- roc_all[, .(avg_mean = mean(auc_asym, na.rm = TRUE),
                         avg_sd   = sd(auc_asym,   na.rm = TRUE),
                         n_pdb    = mean(n_pdb)),
                     by = .(potential, split)]
setorder(aggregati, potential, n_pdb)

# ---- One plot per potential ----
make_plot <- function(dat) {
  ggplot(dat, aes(x = n_pdb, y = avg_mean)) +
    geom_line(color = "#CC4678FF", linewidth = 1) +
    geom_point(color = "#CC4678FF", size = 2) +
    geom_errorbar(aes(ymin = avg_mean - avg_sd, ymax = avg_mean + avg_sd),
                  width = diff(range(dat$n_pdb)) * 0.01,
                  color = "#CC4678FF", alpha = 0.6) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
    labs(
      x = TeX("Number of PDB files used to compute the statistical potentials"),
      y = TeX("Average ROC AUC of the asymmetric potential ($\\pm \\sigma$)")
    ) +
    tema_bianco
}

plot_gr  <- make_plot(aggregati[potential == "g(r)"])
plot_cdr <- make_plot(aggregati[potential == "CDR"])

print(plot_gr)
print(plot_cdr)

ggsave(file.path(base, "roc_auc_asym_convergence_GR.png"),
       plot_gr,  width = 7, height = 5, dpi = 300, bg = "white")
ggsave(file.path(base, "roc_auc_asym_convergence_CDR.png"),
       plot_cdr, width = 7, height = 5, dpi = 300, bg = "white")

roc_all <- rbindlist(list(
  load_roc("convergence_asymmetric_GR",  "g(r)"),
  load_roc("convergence_asymmetric_CDR", "CDR")
))
