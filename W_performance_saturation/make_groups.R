################################################################################
### THIS SHALL BE RUN ONCE, AT THE BEGINNING
################################################################################

library(data.table)

pdb_dir  <- "/Users/lorenzosisti/Downloads/database_settembre_renamed/"
out_root <- "/Users/lorenzosisti/TiNDER/data/convergence_asymmetric_CDR/"

seed_list <- c(101, 123, 202, 303, 404, 505, 606, 707, 808, 909)
partitioning_numbers <- c(2, 3, 5, 8, 13, 21, 34, 55, 89, 144)
ONLY_SELECTED_GROUP  <- TRUE   # FALSE = all n groups; TRUE = 1 random group

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
