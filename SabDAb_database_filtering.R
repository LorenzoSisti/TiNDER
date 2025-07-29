### SCRIPT USED TO OBTAIN A CLEAN, NON REDUNDANT ANTIBODY-ANTIGEN COMPLEXES DATASET ###
### THE SCRIPT RETURNS FILES CONTAINING 2 ANTIBODY CHAINS (HEAVY AND LIGHT) WITH COMPLETE CDRs AND A SINGLE-CHAIN ANTIGEN ###
### THE FILTERING IS DIVIDED IN 3 STEPS: (1) GENERATING FASTAS (2) RUNNING CD-HIT TO AVOID REDUNDANCY (3) CREATING A NON REDUNDANT DATASET ### 

### Libraries ###
library(bio3d)
library(tibble)
library(dplyr)
library(readr)
library(purrr)

### Setting paths and directories ###
work_dir <- "/path/to/your/work/directory/" # Main directory
imgt_dir <- file.path(work_dir, "imgt") # This directory contains the imgt-named structure collection downloaded from SabDAb 
unique_bound_dir <- file.path(work_dir, "unique_bound") # This directory contains unique, bound PDB structured with renamed chains (H, L, A)
dest_dir <- file.path(work_dir, "non_redundant_complete") # Destination directory of all the filtering process

### MAIN 1: GENERATING FASTAS  ###

#From SabDAb, download the .tsv file containing informations on all the antibody-antigen complexes (https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/search/?all=true#downloads)
#By operating on the .tsv, let's discard inappropriate PDBs
df_unique <- read_tsv(file.path(work_dir, "sabdab_summary_all.tsv"), show_col_types = FALSE) %>%
  select(pdb, Hchain, Lchain, antigen_chain, antigen_type, scfv) %>%
  filter(
    !scfv, # Discard the scFv antibodies
    antigen_type %in% c("protein", "peptide"), # Keep just the files whose ligand is a peptide or a protein
    !is.na(Lchain) & !is.na(Hchain), # Discard the files with just the H or the L chain
    nchar(antigen_chain) == 1 # Discard the files with more than one ligand
  ) %>%
  distinct(pdb, .keep_all = TRUE) # If a PDB file has more than a complex in the asymmetric unit due to non-crystallographic symmetry (NCS), keep just one


### Define the filtering function ###

filter_dataset <- function(pdb, Hchain, Lchain, antigen_chain) {
  file_path <- file.path(imgt_dir, paste0(pdb, ".pdb")) # Define the PDB path
  
  if (!file.exists(file_path)) {
    message("File not found: ", file_path)
    return(NULL)
  }
  
  pdb_obj <- read.pdb(file_path) # Read the PDB file
  
  # Select the chains
  ch_h_ndx <- atom.select(pdb_obj, chain = Hchain)
  ch_l_ndx <- atom.select(pdb_obj, chain = Lchain)
  ch_ag_ndx <- atom.select(pdb_obj, chain = antigen_chain)
  
  # Combine the chains
  all_ndx <- combine.select(ch_h_ndx, ch_l_ndx, ch_ag_ndx, operator = 'OR')
  trimmed_pdb <- trim(pdb_obj, all_ndx)
  
  # Rename the heavy and light chains as H and L, which will be useful for Rosetta analysis later on
  trimmed_pdb$atom$chain[trimmed_pdb$atom$chain == Hchain] <- "H"
  trimmed_pdb$atom$chain[trimmed_pdb$atom$chain == Lchain] <- "L"
  trimmed_pdb$atom$chain[trimmed_pdb$atom$chain == antigen_chain] <- "A"
  
  # Save the unique, bound, renamed PDB files
  # The PDBs present in this directory will be further processed and filtered  
  new_filename <- pdb
  output_path <- file.path(unique_bound_dir, paste0(new_filename, ".pdb"))
  write.pdb(trimmed_pdb, output_path)
  
  # Extract the sequences. This is necessary to run the CD-HIT sequence identity and clustering analysis
  seq_H <- pdbseq(pdb_obj, atom.select(pdb_obj, chain = Hchain, "calpha"))
  seq_L <- pdbseq(pdb_obj, atom.select(pdb_obj, chain = Lchain, "calpha"))
  seq_ag <- pdbseq(pdb_obj, atom.select(pdb_obj, chain = antigen_chain, "calpha"))
  
  # Return a tibble with all the sequences
  tibble(
    pdb_id = new_filename,
    ch_h = paste0(unname(seq_H), collapse = ""),
    ch_l = paste0(unname(seq_L), collapse = ""),
    ch_ag = paste0(unname(seq_ag), collapse = "")
  )
}

# Apply the functon to the whole dataset
results <- pmap_dfr( # Applies a function to multiple columns (like a loop) and returns a data frame
  list(df_unique$pdb, df_unique$Hchain, df_unique$Lchain, df_unique$antigen_chain),
  function(pdb, H, L, ag) {
    filter_dataset(pdb = pdb, Hchain = H, Lchain = L, antigen_chain = ag)
  }
)

# Save the sequences in .csv format
write_csv(results, file.path(work_dir, "filtered_sequences.csv"))

# Write the FASTA files which will be given to CD-HIT
# The FASTA will have H chain, L chain and antigen chain concatenated (as they are when using CD-HIT from SabDAb)
fasta_lines <- c()
for (i in 1:nrow(results)) {
  header <- paste0(">", results$pdb_id[i])
  sequence <- paste0(results$ch_h[i], results$ch_l[i], results$ch_ag[i])
  fasta_lines <- c(fasta_lines, header, sequence)
}

writeLines(fasta_lines, file.path(work_dir, "filtered_sequences.fasta"))

### MAIN 2: AVOIDING REDUNDANCY (from bash) ###

### We now have to filter the unique, bound, renamed FASTA files by sequence identity ###
### From bash in linux or MacOS, download CD-HIT as in https://github.com/weizhongli/cdhit/wiki and run:

# cd-hit -i filtered_sequences.fasta -o non_redundant_sequences.fasta 
# grep "^>" non_redundant_sequences.fasta | sed 's/^>//' > df_dataset_clusters.csv

### Clustering is run with default CD-HIT parameters (90% seq. id.) ###

### MAIN 3: CREATING A NON REDUNDANT COLLECTION WITH COMPLETE CDRs  ###

### Load the csv file from CD-HIT with the non-redundant set of proteins
post_cd_hit_pdb <- read.csv(file.path(work_dir, "df_dataset_clusters.csv"))
colnames(post_cd_hit_pdb) <- "pdb" # Rename the PDB column

### Load the tsv file downloaded from SAbDAb with the collection of bound Ab:::Ag complexes that have complete CDRs according to IMGT nomenclature
### https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/cdrsearch/?CDRdef=IMGT&CDRtype=All&completeCDR=Complete&CDRlength=&agtype=All&method=All&species=All&resolution=&rfactor=&ltype=All&CDRseq=#downloads
complete_cdrs <- read_tsv((file.path(work_dir, "complete_cdrs_imgt_summary.tsv")))

pdb_ids <- intersect(post_cd_hit_pdb$pdb, complete_cdrs$pdb)

# For every PDB file
for (pdb_id in pdb_ids) {
  pdb_filename <- paste0(pdb_id, ".pdb")
  source_path <- file.path(unique_bound_dir, pdb_filename)
  destination_path <- file.path(dest_dir, pdb_filename)
  
  # If the file exists, copy it in the destination directory
  if (file.exists(source_path)) {
    file.copy(from = source_path, to = destination_path, overwrite = FALSE)
  } else {
    message("File not found: ", pdb_filename)
  }
}
