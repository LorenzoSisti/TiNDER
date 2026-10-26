### THIS SCRIPT PERFORMS CLUSTERING ANALYSIS ON ANTIBODY-ANTIGEN COMPLEXES INTERFACE DESCRIPTORS. ###
### IT TAKES A PRE-COMPUTED DESCRIPTORS DATAFRAME AS INPUT###
### IT IS DIVIDED INTO THE FOLLOWING MAIN BLOCKS: ###
### 1) DATA IMPORT AND NORMALIZATION ###
### 2) CLUSTERING OF THE INTERFACE DESCRIPTORS ###
### 3) PRINCIPAL COMPONENT ANALYSIS (PCA) ###
### 4) PCA-BASED CLUSTERING AND CLUSTER CENTROID EXTRACTION FOR DOCKING EXPERIMENTS ###

### Libraries ###
library(bio3d)
library(stats)
library(cluster)
library(hopkins)
library(corrr)
library(ggplot2)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(corrplot)
set.seed(1234)

### Function to find the nearest PDB to a cluster centroid
find_cluster_centroid <- function(cl) {
  # Subset only the PDB entries in this cluster
  cluster_data <- reduced_data[clusters == cl, , drop = FALSE]
  cluster_names <- rownames(cluster_data)
  
  # Calculate the centroid
  centroid_coords <- colMeans(cluster_data)
  
  # Compute Euclidean distance of each PDB to the centroid
  dists <- apply(cluster_data, 1, function(row) {
    sqrt(sum((row - centroid_coords)^2))
  })
  
  # Identify the PDB closest to the centroid
  pdb_nearest <- cluster_names[which.min(dists)]
  
  # Return as a small data frame
  data.frame(
    Cluster = cl,
    Centroide_PDB = pdb_nearest,
    Distanza = min(dists)
  )
}

#------------------------------------------
# DATA IMPORT AND NORMALIZATION
#------------------------------------------

# Load CSV with all calculated (non-normalized) descriptors
descriptors <- read.csv("/Users/lorenzosisti/Downloads/interface_descriptors_nuovo/merged_descriptors.csv")

# Convert descriptor data frame to numeric matrix (excluding PDB_ID column)
descriptor_matrix <- as.matrix(descriptors[, - which(names(descriptors) == "PDB_ID")])
rownames(descriptor_matrix) <- descriptors$PDB_ID
scaled_data <- scale(descriptor_matrix) # Normalize each column using Z-score (mean=0, sd=1)

#------------------------------------------
# CLUSTER ANALYSIS ON THE WHOLE DESCRIPTOR DATASET
#------------------------------------------

### Clustering tendency assessment

hopkins(scaled_data, m = ceiling(nrow(scaled_data)/10)) # Statistical method: Hopkins statistic
fviz_dist(dist(scaled_data), show_labels = FALSE) # Visual method: distance heatmap

### Clustering based on euclidean distance between observtions

# Partitioning method: k-means with silhouette evaluation
k_means_partitioning <- fviz_nbclust(scaled_data, kmeans, method = "silhouette", k.max = 50) 

# Hierarchical method (AGNES, average linkage)
hierarchical_partitioning <- agnes(scaled_data, metric = "euclidean", stand = FALSE, method = "average")
fviz_dend(hierarchical_partitioning, show_labels = FALSE)

cophenetic_distance <- cophenetic(hierarchical_partitioning) # Evaluate goodness of clustering: cophenetic correlation
cor(dist(scaled_data, method = "euclidean"), cophenetic_distance)

# In our case, average linkage shows highest cophenetic correlation and was therefore chosen as the default linking method (this is true for all the clustering analysis)

### Clustering based on Pearson correlation distance

# Partitioning method: k-means with Pearson-based distance
k_means_partitioning_pearson <- fviz_nbclust(scaled_data, kmeans, method = "silhouette", diss = get_dist(scaled_data, method = "pearson"), k.max = 50) #La matrice delle distanze tra scaled_data è creata usando la metrica euclidea

# Hierarchical method (average linkage)
# We don't use AGNES anymore because it does not support correlation-based distances
dissimilarity_matrix <- get_dist(scaled_data, method = "pearson")
hierarchical_partitioning_pearson <- hclust(dissimilarity_matrix, method = "average")
fviz_dend(hierarchical_partitioning_pearson, show_labels = FALSE)

# Cophenetic correlation for Pearson-based clustering
cor(get_dist(scaled_data, method = "pearson"), cophenetic(hierarchical_partitioning_pearson))


#------------------------------------------
# BETWEEN-DESCRIPTORS CORRELATION ANALYSIS
#------------------------------------------

library(viridis)

# 1. Define the desired order of descriptors

size_descriptors <- c(
  "antigen_length",
  "interface_radius",
  "paratope_aa_number",
  "epitope_aa_number",
  "interface_gyration_radius",
  "sc_int_area"
)

physicochemical_descriptors <- c(
  "avg_interface_Hr",
  "sc",
  "ec_avg"
)

topological_descriptors <- c(
  "betweenness_centrality",
  "closeness_centrality",
  "shortest_path",
  "degree",
  "clustering_coefficient",
  "density",
  "modularity",
  "graph_energy",
  "graph_energy_of_complexes"
)

descriptor_order <- c(
  size_descriptors,
  physicochemical_descriptors,
  topological_descriptors
)

# 2. Reorder descriptors and calculate the correlation matrix

ordered_descriptors <- descriptor_matrix[
  , descriptor_order, drop = FALSE
]

cor_matrix <- cor(
  ordered_descriptors,
  use = "pairwise.complete.obs",
  method = "pearson"
)

# 3. Create a numbered correspondence table

descriptor_labels <- data.frame(
  Number = seq_along(descriptor_order),
  Category = c(
    rep("Size", length(size_descriptors)),
    rep("Physicochemical",
        length(physicochemical_descriptors)),
    rep("Topological", length(topological_descriptors))
  ),
  Descriptor = descriptor_order
)

print(descriptor_labels, row.names = FALSE)

# 4. Replace axis labels with numbers from 1 to 18

dimnames(cor_matrix) <- list(
  as.character(seq_len(ncol(cor_matrix))),
  as.character(seq_len(ncol(cor_matrix)))
)

# 5. Create a reversed Plasma palette
#    High positive correlations (+1) = purple


# Plasma palette, with purple for high positive correlations
plasma_col <- colorRampPalette(
  viridis::plasma(256, direction = -1)
)(200)

# Increase transparency: alpha = 0.65
plasma_transparent <- grDevices::adjustcolor(
  plasma_col,
  alpha.f = 0.75
)

corrplot(
  cor_matrix,
  method = "color",
  type = "lower",
  order = "original",
  col = plasma_transparent,
  col.lim = c(-1, 1),
  
  # Color scale
  cl.pos = "b",
  cl.length = 5,
  cl.ratio = 0.15,
  cl.cex = 0.8,
  
  # Labels
  tl.pos = "l",
  tl.col = "black",
  tl.cex = 0.9,
  tl.srt = 0,
  
  # Appearance
  addgrid.col = "white",
  diag = TRUE,
  mar = c(1, 5, 3, 1)
)


install.packages(c("igraph", "ggraph", "ggplot2", "dplyr"))

library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)


#------------------------------------------
# CIRCULAR NETWORK OF DESCRIPTOR CORRELATIONS
#------------------------------------------

# 1. Define the correlation threshold
cor_threshold <- 0.50

# 2. Define nodes and categories
nodes <- descriptor_labels %>%
  mutate(
    name = as.character(Number),
    Category = factor(
      Category,
      levels = c("Size", "Physicochemical", "Topological")
    )
  )

# 3. Extract unique pairs of descriptors
#    Exclude the diagonal and retain strong correlations

cor_df <- as.data.frame(as.table(cor_matrix))
colnames(cor_df) <- c("from", "to", "correlation")

edges <- cor_df %>%
  mutate(
    from = as.character(from),
    to = as.character(to)
  ) %>%
  filter(
    as.integer(from) < as.integer(to),
    abs(correlation) >= cor_threshold
  ) %>%
  mutate(
    sign = ifelse(correlation > 0, "Positive", "Negative"),
    strength = abs(correlation)
  )

# 4. Build the undirected network
g <- graph_from_data_frame(
  edges,
  directed = FALSE,
  vertices = nodes
)

# 5. Plot the circular graph
set.seed(1234)

ggraph(g, layout = "circle") +
  geom_edge_link(
    aes(
      edge_colour = sign,
      edge_width = strength
    ),
    alpha = 0.65,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(shape = Category),
    size = 7,
    color = "black",
    fill = "black",
    stroke = 1
  ) +
  scale_shape_manual(
    name = "Descriptor category",
    values = c(
      "Size" = 21,             # Circle
      "Physicochemical" = 22,  # Square
      "Topological" = 24       # Triangle
    ),
    breaks = c("Size", "Physicochemical", "Topological")
  ) +
  geom_node_text(
    aes(label = name),
    size = 3,
    color = "white",
    fontface = "bold"
  ) +
  scale_edge_colour_manual(
    name = "Correlation sign",
    values = c(
      "Positive" = "#7B3294",
      "Negative" = "#E58B55"
    )
  ) +
  scale_edge_width(
    name = "|Correlation|",
    range = c(0.3, 2)
  ) +
  coord_fixed() +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ))


#------------------------------------------
# PCA AND CLUSTER ANALYSIS ON PRINCIPAL COMPONENTS
#------------------------------------------

data_pca <- princomp(scaled_data)
summary(data_pca)

fviz_eig(
  data_pca,
  addlabels = FALSE,
  ncp = 18,
  geom = "bar",
  barfill = "#CC4678FF",
  barcolor = "#CC4678FF"
) +
  labs(title = NULL)
reduced_data <- data_pca$scores[, 1:5] # Keep first 5 principal components

ggplot(as.data.frame(reduced_data),
       aes(x = Comp.1, y = Comp.2)) +
  geom_point(color = "#CC4678", alpha = 0.6, size = 2) +
  labs(
    title = NULL,
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal()

# Re-assess clustering tendency on PCA-reduced data
hopkins(reduced_data, m = ceiling(nrow(scaled_data)/10)) # Statistichal method
fviz_dist(dist(reduced_data), show_labels = FALSE)

### Clustering based on euclidean distance between observtions

# K-means
pca_k_means_partitioning <- fviz_nbclust(reduced_data, kmeans, method = "silhouette", k.max = 50)

# Hierarchical clustering (AGNES, average linkage)
pca_hierarchical_partitioning <- agnes(reduced_data, metric = "euclidean", stand = FALSE, method = "average")
fviz_dend(pca_hierarchical_partitioning, show_labels = FALSE)
pca_cophenetic_distance <- cophenetic(pca_hierarchical_partitioning) # Cophenetic correlation
cor(dist(reduced_data, method = "euclidean"), pca_cophenetic_distance)

### Clustering based on Pearson correlation distance

# K-means
pca_k_means_partitioning_pearson <- fviz_nbclust(reduced_data, kmeans, method = "silhouette", diss = get_dist(reduced_data, method = "pearson"), k.max = 50) #La matrice delle distanze tra scaled_data è creata usando la metrica euclidea

# Hierarchical (we don't use AGNES here because it doesn't support correlation-based distances)
pca_dissimilarity_matrix <- get_dist(reduced_data, method = "pearson")
pca_hierarchical_partitioning_pearson <- hclust(pca_dissimilarity_matrix, method = "average")
fviz_dend(pca_hierarchical_partitioning_pearson, show_labels = FALSE)
cor(get_dist(reduced_data, method = "pearson"), cophenetic(pca_hierarchical_partitioning_pearson))

#------------------------------------------
# EXTRACTING PDB FILES FOR DOCKING ANALYSIS
#------------------------------------------

# Cut the Pearson-based dendrogram into 50 clusters
clusters <- cutree(pca_hierarchical_partitioning_pearson, k = 50)

# Find centroid representatives for each cluster
cluster_centroid <- lapply(unique(clusters), find_cluster_centroid)
cluster_centroid_df <- do.call(rbind, cluster_centroid)


### N.B. FORSE DEVI CORREGGERE UN ERRORE: NELLA FUNZIONE DI ESTRAZIONE SEQUENZE SI ESTRAE DAL DB WHOLE E NON PCA


#------------------------------------------
# ESTRAZIONE SEQUENZE FASTA DAI FILE LOCALI
#------------------------------------------

# 1. CARTELLA INPUT: Dove si trovano i tuoi file PDB processati e tagliati
cartella_pdb_locali <- "/Users/lorenzosisti/Downloads/database_settembre_renamed/"

# 2. CARTELLA OUTPUT: Dove salvare i singoli file FASTA per AlphaFold 3
cartella_destinazione <- "/Users/lorenzosisti/Downloads/FASTAS_FOR_AF3_8_giugno/"

# Crea la cartella di destinazione se non esiste ancora
if (!dir.exists(cartella_destinazione)) {
  dir.create(cartella_destinazione, recursive = TRUE)
}

# Estrazione dei PDB centroidi unici identificati dal tuo clustering
pdb_selezionati <- unique(cluster_centroid_df$Centroide_PDB)
message("Inizio estrazione sequenze. Saranno creati ", length(pdb_selezionati), " file FASTA individuali...")

for (pdb_id in pdb_selezionati) {
  
  # Pulizia del nome da spazi bianchi
  pdb_clean <- trimws(pdb_id)
  
  # Gestione Case-Sensitivity per il file PDB di input
  nome_file_pdb <- paste0(tolower(pdb_clean), ".pdb") 
  percorso_completo_pdb <- file.path(cartella_pdb_locali, nome_file_pdb)
  
  if (!file.exists(percorso_completo_pdb)) {
    nome_file_pdb <- paste0(pdb_clean, ".pdb")
    percorso_completo_pdb <- file.path(cartella_pdb_locali, nome_file_pdb)
    
    if (!file.exists(percorso_completo_pdb)) {
      warning(paste("File PDB NON trovato nella cartella locale:", percorso_completo_pdb))
      next # Salta al PDB successivo se questo manca
    }
  }
  
  # Lettura del PDB locale
  pdb_obj <- tryCatch({
    read.pdb(percorso_completo_pdb, verbose = FALSE)
  }, error = function(e) {
    warning(paste("Errore nella lettura del file:", percorso_completo_pdb))
    return(NULL)
  })
  
  # Se il file è valido, procediamo alla creazione del singolo FASTA
  if (!is.null(pdb_obj)) {
    
    # --- NOVITÀ: Definiamo il percorso del singolo file FASTA (es. 7x8t.fasta) ---
    nome_file_fasta <- paste0(tolower(pdb_clean), ".fasta")
    percorso_completo_fasta <- file.path(cartella_destinazione, nome_file_fasta)
    
    # Inizializza (o svuota se già esistente) il file FASTA specifico per questo PDB
    cat("", file = percorso_completo_fasta, append = FALSE)
    
    message("Scrivo il file FASTA: ", nome_file_fasta)
    
    # Estraiamo le catene del PDB
    catene <- unique(pdb_obj$atom$chain)
    
    for (catena in catene) {
      if (catena == "" || is.na(catena)) next
      
      # Isola la catena ed estrae la sequenza
      pdb_isolato <- trim.pdb(pdb_obj, chain = catena)
      sequenza_vettore <- pdbseq(pdb_isolato)
      
      if (length(sequenza_vettore) > 0) {
        # Converte il vettore in stringa continua
        sequenza_stringa <- paste(sequenza_vettore, collapse = "")
        
        # Header FASTA (>pdb_id|Chain_X)
        header <- paste0(">", tolower(pdb_clean), "|Chain_", catena)
        
        # Scrive la catena nel file FASTA specifico di questo PDB
        cat(header, "\n", sequenza_stringa, "\n", file = percorso_completo_fasta, append = TRUE)
      }
    }
  }
}

message("Processo completato! Trovi i 50 file .fasta pronti dentro: ", cartella_destinazione)
