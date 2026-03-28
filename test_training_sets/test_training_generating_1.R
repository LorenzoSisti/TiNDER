### Required libraries
library(bio3d)
library(dplyr)
library(tidyr)
library(future)
library(furrr)

# Define the path to a custom function files
#source("/path/to/your/custom/functions/files/functions.R")
source("/Users/lorenzosisti/Documents/Script_ottimizzati_funzioni/functions.R")


### Set up parallelization to speed up computation
plan(multisession, workers = parallel::detectCores() - 1)

### Define directories and global parameters
pdb_dir <- "/Users/lorenzosisti/Downloads/database_settembre_renamed/"
docked_dir <- "/Users/lorenzosisti/Downloads/models/"
results_dir <- "/Users/lorenzosisti/Downloads/potenziali_statistici_test_training"
dir.create(results_dir, showWarnings = FALSE)

# Distance cutoff (Å) to define contact between side-chains centroids
DistCutoff <- 8.5
Nsteps <- 3

# Amino acids ordered according to the Kyte-Doolittle hydrophobicity scale
amino_acids <- c("ARG", "LYS", "ASN", "ASP", "GLN", "GLU", "HIS", "PRO", "TYR", "TRP",
                 "SER", "THR", "GLY", "ALA", "MET", "CYS", "PHE", "LEU", "VAL", "ILE")

# Retrieve all PDB files from the docking folder
all_pdbs <- list.files(pdb_dir, pattern = "*.pdb", recursive = TRUE, full.names = TRUE)
all_docked_pdbs <- list.files(docked_dir, pattern = "*.pdb", recursive = TRUE, full.names = TRUE)


# Imposto il seed per la riproducibilità
set.seed(123) 

# Estraggo solo il nome del file, rimuovendo il percorso completo e l'estensione ".pdb"
# Es: "/.../9lbsa.pdb" diventa "9lbsa"
nomi_base <- sub("\\.pdb$", "", basename(all_pdbs))

# Trovo la posizione (indici) di tutti i file che hanno esattamente 4 caratteri nel nome base
indici_4_char <- which(nchar(nomi_base) == 4)

# Estraggo casualmente 187 indici ESCLUSIVAMENTE tra quelli che hanno 4 caratteri
indici_187 <- sample(indici_4_char, 187)
  
# Creo i due sottogruppi
pdbs_187  <- all_pdbs[indici_187]   # Le 187 proteine scelte a caso (solo 4 caratteri)
pdbs_2000 <- all_pdbs[-indici_187]  # Le rimanenti 2000 proteine

# 1. Estraggo i codici PDB (i 4 caratteri) dalla lista pdbs_187
codici_187 <- sub("\\.pdb$", "", basename(pdbs_187))

# 2. Filtro all_docked_pdbs
# Prendo solo i file in cui i primi 4 caratteri del nome file corrispondono a uno dei codici_187
docked_187 <- all_docked_pdbs[substr(basename(all_docked_pdbs), 1, 4) %in% codici_187]

# Definisco i percorsi delle due nuove sottocartelle
test_dir <- file.path(results_dir, "test_dir")
training_dir <- file.path(results_dir, "training_dir")

# Creo le sottocartelle (showWarnings = FALSE evita errori se esistono già)
dir.create(test_dir, showWarnings = FALSE)
dir.create(training_dir, showWarnings = FALSE)

# Copio i file di test (docked_187) nella test_dir
# overwrite = TRUE sovrascrive i file se dovessi far girare lo script più volte
file.copy(from = docked_187, to = test_dir, overwrite = TRUE)

# Copio i file di training (pdbs_2000) nella training_dir
file.copy(from = pdbs_2000, to = training_dir, overwrite = TRUE)
