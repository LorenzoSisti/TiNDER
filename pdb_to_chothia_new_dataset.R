# ============================================================
# Renumber H/L chains of docking poses using Chothia reference
# models/ contains pdbid.pdb (chains H and L, Chothia numbering)
# ============================================================
library(bio3d)
library(dplyr)
library(tibble)
library(purrr)
library(stringr)

# ---------- Parameters ----------
poses_dir     <- "/Users/lorenzosisti/Downloads/models_2_renamed"
ref_dir       <- "/Users/lorenzosisti/Downloads/database_settembre_renamed"
out_dir       <- "/Users/lorenzosisti/Downloads/models_2_renamed_chothia"
chains        <- c("H", "L")
keep_unmapped <- FALSE   # TRUE: keep residues outside the overlap with old numbers
# FALSE: drop them

# ---------- Helpers ----------

# One row per residue (file order) for a given chain
residue_table <- function(pdb, chain_id) {
  at  <- pdb$atom
  idx <- which(at$chain == chain_id)
  if (length(idx) == 0) return(NULL)
  at  <- at[idx, ]
  ins <- ifelse(is.na(at$insert), "", at$insert)
  key <- paste(at$chain, at$resno, ins, sep = "|")
  first <- !duplicated(key)
  tibble(
    key     = key[first],
    resno   = at$resno[first],
    insert  = ins[first],
    resname = at$resid[first]
  ) %>%
    mutate(aa = vapply(resname, function(x) {
      a <- tryCatch(aa321(x), error = function(e) "X")
      if (length(a) == 0 || is.na(a)) "X" else unname(a)
    }, character(1)))
}

seq_of <- function(tab) paste(tab$aa, collapse = "")

# Exact (zero tolerance) location of one sequence inside the other
locate <- function(pose_seq, ref_seq) {
  pos <- str_locate_all(pose_seq, fixed(ref_seq))[[1]]
  if (nrow(pos) == 1) return(list(type = "ref_in_pose", start = pos[1, "start"]))
  if (nrow(pos) > 1)  return(list(type = "ambiguous",   start = NA))
  pos <- str_locate_all(ref_seq, fixed(pose_seq))[[1]]
  if (nrow(pos) == 1) return(list(type = "pose_in_ref", start = pos[1, "start"]))
  if (nrow(pos) > 1)  return(list(type = "ambiguous",   start = NA))
  list(type = "no_match", start = NA)
}

# Mapping pose residue key -> new resno / insert (NULL if no exact match)
map_chain <- function(pose_tab, ref_tab) {
  loc <- locate(seq_of(pose_tab), seq_of(ref_tab))
  if (!loc$type %in% c("ref_in_pose", "pose_in_ref"))
    return(list(map = NULL, status = loc$type))
  
  n_pose <- nrow(pose_tab); n_ref <- nrow(ref_tab); s <- loc$start
  if (loc$type == "ref_in_pose") {
    pose_i <- s:(s + n_ref - 1);  ref_i <- seq_len(n_ref)
  } else {
    pose_i <- seq_len(n_pose);    ref_i <- s:(s + n_pose - 1)
  }
  list(
    map = tibble(pose_key   = pose_tab$key[pose_i],
                 new_resno  = ref_tab$resno[ref_i],
                 new_insert = ref_tab$insert[ref_i]),
    status = loc$type
  )
}

# ---------- Indexes ----------
ref_index <- tibble(path = list.files(ref_dir, pattern = "\\.pdb$",
                                      full.names = TRUE, ignore.case = TRUE)) %>%
  mutate(pdbid = tolower(str_remove(basename(path), "\\.pdb$")))

if (anyDuplicated(ref_index$pdbid))
  stop("Duplicated pdb ids in ", ref_dir)

pose_index <- tibble(path = list.files(poses_dir, pattern = "\\.pdb$",
                                       full.names = TRUE, recursive = TRUE,
                                       ignore.case = TRUE)) %>%
  mutate(
    file      = basename(path),
    algorithm = basename(dirname(path)),
    pdbid     = tolower(sub("_.*$", "", file))
  )

cat("Poses found:     ", nrow(pose_index), "\n")
cat("References found:", nrow(ref_index), "\n")

# ---------- Core function ----------
renumber_pose <- function(pose_path, algorithm, pdbid) {
  
  res <- list(file = basename(pose_path), algorithm = algorithm, pdbid = pdbid,
              status = NA_character_, detail = NA_character_)
  fail <- function(msg) { res$status <- "FAILED"; res$detail <- msg; res }
  
  ref_path <- ref_index$path[ref_index$pdbid == pdbid]
  if (length(ref_path) == 0) return(fail("no reference with this pdb id"))
  
  pose <- tryCatch(read.pdb(pose_path, rm.alt = TRUE, verbose = FALSE),
                   error = function(e) NULL)
  ref  <- tryCatch(read.pdb(ref_path,  rm.alt = TRUE, verbose = FALSE),
                   error = function(e) NULL)
  if (is.null(pose)) return(fail("cannot read pose"))
  if (is.null(ref))  return(fail("cannot read reference"))
  
  # Map each chain; reject the whole file if any chain fails
  maps <- list()
  for (ch in chains) {
    pt <- residue_table(pose, ch)
    rt <- residue_table(ref,  ch)
    if (is.null(pt)) return(fail(paste("pose lacks chain", ch)))
    if (is.null(rt)) return(fail(paste("reference lacks chain", ch)))
    m <- map_chain(pt, rt)
    if (is.null(m$map)) return(fail(paste0("chain ", ch, ": ", m$status)))
    maps[[ch]] <- m$map
  }
  all_map <- bind_rows(maps)
  
  # Apply mapping at atom level
  at   <- pose$atom
  ins  <- ifelse(is.na(at$insert), "", at$insert)
  akey <- paste(at$chain, at$resno, ins, sep = "|")
  
  hit <- match(akey, all_map$pose_key)
  ok  <- !is.na(hit)
  
  new_resno  <- at$resno
  new_insert <- ins
  new_resno[ok]  <- all_map$new_resno[hit[ok]]
  new_insert[ok] <- all_map$new_insert[hit[ok]]
  
  keep <- rep(TRUE, nrow(at))
  if (!keep_unmapped) keep <- !(at$chain %in% chains) | ok
  
  at$resno  <- new_resno
  at$insert <- ifelse(new_insert == "", NA, new_insert)
  
  # Subset atoms and matching xyz triplets
  keep_idx <- which(keep)
  xyz_idx  <- as.vector(rbind(3 * keep_idx - 2, 3 * keep_idx - 1, 3 * keep_idx))
  
  out_pdb        <- pose
  out_pdb$atom   <- at[keep_idx, ]
  out_pdb$xyz    <- pose$xyz[, xyz_idx, drop = FALSE]
  class(out_pdb$xyz) <- "xyz"
  out_pdb$calpha <- out_pdb$atom$elety == "CA"
  out_pdb$atom$eleno <- seq_len(nrow(out_pdb$atom))
  
  out_path <- file.path(out_dir, algorithm, basename(pose_path))
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  write.pdb(out_pdb, file = out_path)
  
  res$status <- "OK"
  res$detail <- sprintf("H: %d res, L: %d res renumbered; %d atoms dropped",
                        nrow(maps$H), nrow(maps$L), sum(!keep))
  res
}

# ---------- Run ----------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_df <- pmap_dfr(
  list(pose_index$path, pose_index$algorithm, pose_index$pdbid),
  function(p, a, i) as_tibble(renumber_pose(p, a, i))
)

write.csv(log_df, file.path(out_dir, "renumbering_log.csv"), row.names = FALSE)

# ---------- Summary ----------
log_df %>% count(algorithm, status) %>% print(n = Inf)
log_df %>% filter(status != "OK")   %>% print(n = Inf)
