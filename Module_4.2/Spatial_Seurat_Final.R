##############################################################
# Spatial Transcriptomics – Breast Cancer Visium (CID44971)
# 3rd Statistical Genomics Workshop
# Based on dataset from Wu et al., Nature Genetics 2021
# DOI: 10.1038/s41588-021-00911-1
##############################################################
# Install R packages (only need to do this once, or if you need to update packages)
#install.packages("devtools")
#devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
#install.packages("Seurat")
#install.packages("ggplot2")
#install.packages("SeuratObject", type = "source")
#install.packages("hdf5r")
#install.packages("ape")
#install.packages('BiocManager')
#BiocManager::install('glmGamPoi')
#install.packages("pheatmap")

# >>> INTRO QUESTION:
# Why is spatial transcriptomics important for understanding breast cancer biology?
# Think about tumor heterogeneity, spatial niches, immune infiltration.

##############################################################
# Load required libraries
##############################################################
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ape)   # Moran’s I for spatial autocorrelation
library(pheatmap)
##############################################################
# 1. Load Visium dataset
##############################################################

# Folder containing:
#   filtered_feature_bc_matrix.h5
#   spatial/
data_dir <- "Module_4.2/Data/spatial/"

breast <- Load10X_Spatial(
  data.dir = data_dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay    = "Spatial",
  slice    = "breast"
)

print(breast)
head(breast@meta.data)

# >>> QUESTION 1:
# How many spots and how many genes are present in the dataset?
# Check the "Dimensions" line in breast.

##############################################################
# 2. Quality Control (QC)
##############################################################

# Compute % mitochondrial reads & capture efficiency
breast[["percent.mt"]] <- PercentageFeatureSet(breast, pattern = "^MT-")
breast$cap_eff <- breast$nFeature_Spatial / breast$nCount_Spatial

head(breast@meta.data[, c("nCount_Spatial","nFeature_Spatial","percent.mt","cap_eff")])

# Violin plots of QC metrics
VlnPlot(
  breast,
  features = c("nCount_Spatial","nFeature_Spatial","percent.mt","cap_eff"),
  pt.size = 0.1,
  ncol = 2
)

# Spatial QC visualization
SpatialFeaturePlot(breast, features = "nCount_Spatial")
SpatialFeaturePlot(breast, features = "nFeature_Spatial")
SpatialFeaturePlot(breast, features = "percent.mt")
SpatialFeaturePlot(breast, features = "cap_eff")

# >>> QUESTION 2:
# Do low-quality spots cluster spatially?
# What biological or technical reasons might explain this pattern?

##############################################################
# 3. (Optional) Filter low-quality spots
##############################################################

nFeature_thresh <- quantile(breast$nFeature_Spatial, 0.25)
nCount_thresh   <- quantile(breast$nCount_Spatial, 0.25)

breast <- subset(
  breast,
  subset = nCount_Spatial > nCount_thresh &
    nFeature_Spatial > nFeature_thresh
)

cat("Spots retained:", nrow(breast@meta.data), "\n")

# >>> QUESTION 3:
# After filtering, how many spots remain?
# Based on QC plots — does the filtering seem justified?

##############################################################
# 4. Normalize with SCTransform
##############################################################

breast <- SCTransform(breast, assay = "Spatial", verbose = FALSE)

##############################################################
# 5. Explore marker gene expression
##############################################################

SpatialFeaturePlot(
  breast,
  features = c("ERBB2", "MUC1", "KRT8"),
  pt.size.factor = 2.4
)

# Save high-quality image
plot_krt8 <- SpatialFeaturePlot(
  breast,
  features = "KRT8",
  pt.size.factor = 2.4
) + theme(
  legend.text  = element_text(size = 8),
  legend.title = element_text(size = 12),
  legend.key.size = unit(0.5, "cm")
)

jpeg("Module_4.2/Results/spatial_KRT8_example.jpg", width = 1200, height = 700, quality = 80)
print(plot_krt8)
dev.off()

# >>> QUESTION 4:
# Where are luminal epithelial markers (KRT8, MUC1) expressed?
# Does the pattern correspond to tumor-rich areas?

##############################################################
# 6. Dimensionality Reduction & Clustering
##############################################################

breast <- RunPCA(breast, assay = "SCT")
breast <- FindNeighbors(breast, reduction = "pca", dims = 1:30)
breast <- FindClusters(breast, resolution = 2)
breast <- RunUMAP(breast, dims = 1:30)

# Visualize clusters
p1 <- DimPlot(breast, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(breast, label = TRUE, label.size = 3, pt.size.factor = 2.7)
p1 + p2

# Highlight specific clusters
SpatialDimPlot(
  breast,
  cells.highlight = CellsByIdentities(breast, idents = c(2,1,4,3,5,8)),
  facet.highlight = TRUE,
  ncol = 3,
  pt.size.factor = 3
)

# >>> QUESTION 5:
# How many clusters did you obtain?
# Which clusters look tumor-like? Which appear immune or stromal?

##############################################################
# 7. Differential Expression (Example: Cluster 5 vs 6)
##############################################################

de_markers <- FindMarkers(breast, ident.1 = 5, ident.2 = 6)
head(de_markers)

SpatialFeaturePlot(
  breast,
  features = rownames(de_markers)[1:3],
  ncol = 3,
  pt.size.factor = 2.7,
  alpha = c(0.1,1)
)

# >>> QUESTION 6:
# Do the DE genes reflect meaningful biological differences?
# (e.g., immune vs tumor, ERBB2-high vs ERBB2-low)

##############################################################
# 8. Marker Genes for All Clusters
##############################################################

all_markers <- FindAllMarkers(
  breast,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

head(all_markers)

# >>> QUESTION 7:
# What are the parameters in FindAllMarkers()?
# What statistical test does Seurat use?
# How many positive markers per cluster were identified?

##############################################################
# Heatmap of top 50 markers per cluster
##############################################################

# Top 50 markers per cluster
top50 <- all_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 50)

# Extract SCT scaled data (Seurat v5 syntax)
sct_mat <- GetAssayData(breast, assay = "SCT", layer = "scale.data")

# Select genes present in SCT matrix
heat_genes <- unique(top50$gene)
heat_genes <- intersect(heat_genes, rownames(sct_mat))

# Build matrix for heatmap
heat_mat <- sct_mat[heat_genes, ]

# Order spots by cluster identity
heat_mat <- heat_mat[, order(Idents(breast))]

# Plot heatmap
pheatmap(
  heat_mat,
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_cols = FALSE,
  main = "Top 50 Markers per Cluster"
)


##############################################################
# 9. Spatially Variable Genes (Moran’s I)
##############################################################

breast <- FindSpatiallyVariableFeatures(
  breast,
  selection.method = "moransi",
  assay = "SCT",
  features = VariableFeatures(breast)[1:1000]
)

top.features <- head(SpatiallyVariableFeatures(breast), 6)
top.features

SpatialFeaturePlot(
  breast,
  features = top.features,
  ncol = 3,
  pt.size.factor = 3
)

# >>> FINAL QUESTION:
# How could spatially variable genes help identify:
#   - tumor subclones?
#   - immune infiltration niches?
#   - drug-resistant regions?
##############################################################
# END OF SCRIPT
##############################################################
