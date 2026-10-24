library(data.table)
library(pROC)

### Header: the only part to change
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript roc_batch.R <root_dir> <csv_name> [out_name]")
root     <- args[1]
csv_name <- args[2]
out_name <- if (length(args) >= 3) args[3] else paste0("roc_auc_", basename(root), ".csv")

DOCKQ_THR <- 0.81
DIRECTION <- "<"

### Find all result files
files <- list.files(root, pattern = paste0("^", csv_name, "$"), recursive = TRUE, full.names = TRUE)
message(length(files), " files found")

get_auc <- function(df, col) {
  d <- df[!is.na(df[[col]]) & !is.na(df$DockQ), ]
  if (nrow(d) == 0 || length(unique(d$true_class)) < 2) return(NA_real_)
  as.numeric(auc(roc(d$true_class, d[[col]], direction = DIRECTION, quiet = TRUE)))
}

res <- rbindlist(lapply(files, function(f) {
  df <- as.data.frame(fread(f))
  df$true_class <- ifelse(df$DockQ <= DOCKQ_THR, 1, 0)
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

setorder(res, split, seed, group)
fwrite(res, file.path(root, out_name))
print(res[, .(mean_auc_sym = mean(auc_sym, na.rm = TRUE),
              mean_auc_asym = mean(auc_asym, na.rm = TRUE),
              n_groups = .N), by = split])