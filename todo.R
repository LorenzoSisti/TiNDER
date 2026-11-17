library(bio3d)

docking_models_dir <- "/Users/lorenzosisti/Downloads/models_2"
output_dir         <- "/Users/lorenzosisti/Downloads/models_2_renamed"

chain_metadata <- read.csv("/Users/lorenzosisti/Downloads/metadata_200Ab-Ag_complex.csv")
trimer_models  <- chain_metadata[!chain_metadata$chainL == "", ]   # antibodies only (H + L + A)

all_models <- list.files(docking_models_dir, pattern = "\\.pdb$",
                         recursive = TRUE, full.names = TRUE)

# Keep only the models that belong to an antibody complex
is_antibody <- sapply(basename(all_models), function(f) {
  any(startsWith(f, trimer_models$pdb_wchains))
})
all_models <- all_models[is_antibody]
message(length(all_models), " antibody models to process")

# Sequence of each chain in a bio3d pdb object (named by chain ID)
get_chain_seqs <- function(pdb) {
  sapply(unique(pdb$atom$chain), function(ch) {
    paste(pdbseq(trim.pdb(pdb, chain = ch)), collapse = "")
  })
}

rename_model_chains <- function(pdb_path) {
  
  file_name <- basename(pdb_path)
  
  tryCatch({
    
    # Find the metadata row corresponding to this model
    row <- which(startsWith(file_name, trimer_models$pdb_wchains))
    if (length(row) != 1) stop("no unique metadata match")
    meta <- trimer_models[row, ]
    
    pdb <- read.pdb(pdb_path)
    model_seqs <- get_chain_seqs(pdb)
    if (length(model_seqs) != 3) stop("expected 3 chains, found ", length(model_seqs))
    
    # Exact match for the antibody chains only
    ab_seqs <- c(H = meta$chainH, L = meta$chainL)
    idx <- match(ab_seqs, model_seqs)
    if (anyNA(idx)) stop("no exact sequence match for: ",
                         paste(names(ab_seqs)[is.na(idx)], collapse = ", "))
    
    # The remaining chain is the antigen
    ab_ids <- names(model_seqs)[idx]
    ag_id  <- setdiff(names(model_seqs), ab_ids)
    if (length(ag_id) != 1) stop("could not isolate a single antigen chain")
    
    # Old chain ID -> new chain ID
    chain_map <- setNames(c("H", "L", "A"), c(ab_ids, ag_id))
    pdb$atom$chain <- unname(chain_map[pdb$atom$chain])
    
    # Save, keeping the subfolder structure
    out_path <- file.path(output_dir, sub(docking_models_dir, "", dirname(pdb_path)))
    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
    write.pdb(pdb, file = file.path(out_path, file_name))
    
    # Record whether the antigen was also identical to the metadata
    data.frame(pdb = file_name, ok = TRUE, ag_exact = model_seqs[[ag_id]] == meta$chainA,
               error = NA)
    
  }, error = function(e) {
    message("Error in ", file_name, ": ", e$message)
    data.frame(pdb = file_name, ok = FALSE, ag_exact = NA, error = e$message)
  })
}

log_df <- do.call(rbind, lapply(all_models, rename_model_chains))
write.csv(log_df, file.path(output_dir, "renaming_log.csv"), row.names = FALSE)

table(log_df$ok)
table(log_df$ag_exact, useNA = "ifany")
unique(substr(log_df$pdb[which(log_df$ag_exact == FALSE)], 1, 4))   # PDB IDs with antigen mismatch