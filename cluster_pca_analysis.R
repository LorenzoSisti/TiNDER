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
set.seed(1234)

### Function to find the nearest PDB to a cluster centroid
find_cluster_centroid <- function(cl) {
  # Subset only the PDB entries in this cluster
  cluster_data <- scaled_data[clusters == cl, , drop = FALSE]
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
descriptors <- read.csv("/path/to/your/descriptor/dataframe")

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

# Compute correlation matrix between descriptors
cor_matrix <- cor(descriptor_matrix, use = "pairwise.complete.obs")
print(cor_matrix)

my_col <- colorRampPalette(c("gold1", "white", "dodgerblue2"))(200) # Custom color palette: gold1 → white → dodgerblue2
corrplot(cor_matrix, col = my_col, tl.pos = "n") 
mtext("Interface descriptors correlation matrix", side = 3, line = 3, col = "black", cex = 1.3, font = 2) # Set title

#------------------------------------------
# PCA AND CLUSTER ANALYSIS ON PRINCIPAL COMPONENTS
#------------------------------------------

data_pca <- princomp(scaled_data)
summary(data_pca)

fviz_eig(data_pca, addlabels = FALSE, ncp = 18, geom = "bar") # Scree plot (eigenvalues)
reduced_data <- data_pca$scores[, 1:5] # Keep first 5 principal components

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

# Cut the Pearson-based dendrogram into 25 clusters
clusters <- cutree(pca_hierarchical_partitioning_pearson, k = 25)

# Find centroid representatives for each cluster
cluster_centroid <- lapply(unique(clusters), find_cluster_centroid)
cluster_centroid_df <- do.call(rbind, cluster_centroid)