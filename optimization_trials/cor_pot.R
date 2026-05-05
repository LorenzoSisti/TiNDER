#########################################################################################
# THIS SCRIPT COMPUTES AMINO ACID CONTACT MATRICES AND WHOLE-INTERFACE STATISTICAL POTENTIALS FROM A DIRECTORY OF PDB STRUCTURES. 
# IT ANALYZES ANTIBODY–ANTIGEN INTERACTIONS BY:
# 1) CALCULATING SIDE-CHAIN CENTROIDS AND CONTACT MATRICES
# 2) IDENTIFYING INTERFACE RESIDUES BETWEEN ANTIBODY AND ANTIGEN
# 3) SUMMARIZING CONTACT FREQUENCIES ACROSS ALL STRUCTURES
# 4) COMPUTING ASYMMETRIC AND SYMMETRIC WHOLE-INTERFACE STATISTICAL POTENTIALS
# 5) VISUALIZING THE RESULTS AS HEATMAPS.
#########################################################################################

# IMPLEMENTAZIONI FUTURE: 
# IMPLEMENTARE EQUAZIONE PER SPARSE DATA E VEDERE SE I DUE POTENZIALI (SPARSE VS NON SPARSE) CORRELANO SUFFICIENTEMENTE

### Required libraries
library(bio3d)
library(dplyr)
library(future)
library(furrr)
library(purrr)
library(progressr)
library(pheatmap)
library(patchwork)
library(ggplotify)
library(reshape2)
library(tidyr)

# Define the path to a custom function files
source("/Users/lorenzosisti/TiNDER/functions.R")

### Set up parallelization to speed up computation
plan(multisession, workers = parallel::detectCores() - 1)
handlers(global = TRUE)
handlers("rstudio")

### Define directories and global parameters
pdb_dir <- "/Users/lorenzosisti/Downloads/database_settembre_renamed/"
results_dir <- "/Users/lorenzosisti/TiNDER/optimization_trials/cor_whole_interface_statistical_potential/" 
dir.create(results_dir, showWarnings = FALSE)

non_sparse_asym <- read.csv("/Users/lorenzosisti/TiNDER/optimization_trials/whole_interface_statistical_potential/V_asym_df.csv")
non_sparse_sym <- read.csv("/Users/lorenzosisti/TiNDER/optimization_trials/whole_interface_statistical_potential/V_sym_df.csv")

sparse_asym <- read.csv("/Users/lorenzosisti/TiNDER/optimization_trials/sparse_whole_interface_statistical_potential/V_asym_df.csv")
sparse_sym <- read.csv("/Users/lorenzosisti/TiNDER/optimization_trials/sparse_whole_interface_statistical_potential/V_sym_df.csv")

cor(sparse_asym$value, non_sparse_asym$value)
cor(sparse_sym$value, non_sparse_sym$value)

