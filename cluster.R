# Run the standard workflow for visualization and clustering
seu.combined <- ScaleData(seu.combined, verbose = FALSE)
seu.combined <- RunPCA(seu.combined, npcs = 30, verbose = FALSE)
seu.combined <- RunUMAP(seu.combined, reduction = "pca", dims = 1:15)
seu.combined <- FindNeighbors(seu.combined, reduction = "pca", dims = 1:15)
seu.combined <- FindClusters(seu.combined, resolution = 0.2)

# PCA
pca <- DimPlot(seu.combined, reduction = "pca", group.by = "donor_id") + 
  ggtitle("PCA Plot")
ggsave("figures/PCA.png", pca)

cluster_df <- data.frame(cbind(Idents(seu.combined), metadata$donor_id))
names(cluster_df) <- c("assigned_cluster", "donor_id")
cluster_df$assigned_cluster <- as.numeric(cluster_df$assigned_cluster)
table(cluster_df$assigned_cluster, cluster_df$donor_id)

# Create UMAP
umap_celltype <- DimPlot(seu.combined, reduction = "umap", group.by = "cell_type__custom") + 
  ggtitle("UMAP Plot")
ggsave("figures/UMAP_celltype.png", umap_celltype)

umap_donor <- DimPlot(seu.combined, reduction = "umap", group.by = "donor_id") + 
  ggtitle("UMAP Plot")
ggsave("figures/UMAP_donor.png", umap_donor)

umap_disease <- DimPlot(seu.combined, reduction = "umap", group.by = "disease_status") + 
  ggtitle("UMAP Plot")
ggsave("figures/UMAP_disease.png", umap_disease)

ElbowPlot(object=seu.combined, ndims = 30)
# elbow appears around 15 PCs hence why I used 15 for the umaps 
FeaturePlot(seu.combined, features = c("Braf","Raf1","Mapk1"))
