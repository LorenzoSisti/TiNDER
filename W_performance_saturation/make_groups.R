################################################################################
### THIS SHALL BE RUN ONCE, AT THE BEGINNING
################################################################################

# DECOMMENT THIS IF RUNNING BATCH JOB WITH run_all.sh

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript make_groups.R <rep_root> <rep_id>")
rep_root <- args[1]
rep_id   <- as.integer(args[2])

pdb_dir  <- "/Users/lorenzosisti/Downloads/database_settembre_renamed/"
out_root <- file.path(rep_root, "convergence_asymmetric_GR")   # the CDR path is derived from this

seed_list            <- rep_id * 1000 + 1:10
partitioning_numbers <- c(2, 3, 5, 8, 13, 21, 34, 55, 89, 144)
ONLY_SELECTED_GROUP  <- TRUE

################################################################################
################################################################################
################################################################################

# DECOMMENT THIS IF RUNNING INTERACTIVE JOB 

# library(data.table)
# 
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) < 2) stop("Usage: Rscript make_groups.R <rep_root> <rep_id>")
# rep_root <- args[1]
# rep_id   <- as.integer(args[2])
# 
# pdb_dir  <- "/Users/lorenzosisti/Downloads/database_settembre_renamed/"
# out_root <- file.path(rep_root, "convergence_asymmetric_GR")   # the CDR path is derived from this
# 
# seed_list            <- rep_id * 1000 + 1:10
# partitioning_numbers <- c(2, 3, 5, 8, 13, 21, 34, 55, 89, 144)
# ONLY_SELECTED_GROUP  <- TRUE

pdbL <- sort(list.files(pdb_dir, pattern = "\\.pdb$", recursive = TRUE, full.names = TRUE))

dir.create(file.path(out_root, "groups"), recursive = TRUE, showWarnings = FALSE)
jobs <- list()

for (s in seed_list) for (n in partitioning_numbers) {
  set.seed(s)
  shuf <- sample(pdbL)
  fpg  <- floor(length(shuf) / n); stopifnot(fpg >= 1)
  grp  <- split(shuf[1:(fpg * n)], rep(1:n, each = fpg))
  sel  <- if (ONLY_SELECTED_GROUP) sample(seq_len(n), 1) else seq_len(n)
  
  for (g in sel) {
    list_file <- file.path(out_root, "groups", sprintf("seed_%s_split_%s_group_%s.txt", s, n, g))
    out_dir   <- file.path(out_root, paste0("seed_", s), paste0("split_", n), paste0("group_", g))
    writeLines(grp[[g]], list_file)          # one full PDB path per line
    jobs[[length(jobs) + 1]] <- data.table(seed = s, split = n, group = g,
                                           list_file = list_file, out_dir = out_dir)
  }
}
fwrite(rbindlist(jobs), file.path(out_root, "jobs.tsv"), sep = "\t")
