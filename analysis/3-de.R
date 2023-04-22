# Set levels of seurat object;
# Goal is to compare gene expression between kd and control mice within individual cells
Idents(seu.combined) <- "disease_status"
levels(seu.combined)

# List of cell types
unique(seu.combined$cell_type__custom)

#### Differential Expression Analysis within cell types ####
# Podocyte
marker_genes <- FindMarkers(subset(seu.combined, subset = cell_type__custom == "Podocyte"), 
                             ident.1 = "CTRL", ident.2 = "KDKD")
head(marker_genes %>% arrange(desc(avg_log2FC)), 15)
head(marker_genes, 15)

for (celltype in unique(seu.combined$cell_type__custom)) {
  marker_genes <- FindMarkers(subset(seu.combined, subset = cell_type__custom == celltype), 
                              ident.1 = "CTRL1", ident.2 = "KDKD1")
  print(head(marker_genes, 15))
}