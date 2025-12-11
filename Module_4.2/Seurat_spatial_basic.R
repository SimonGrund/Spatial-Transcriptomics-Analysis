##############################################################
# SPATIAL TRANSCRIPTOMICS WORKSHOP – SEURAT HANDS-ON MODULE
# Dataset: Breast Cancer Visium (CID44971) from the paper: DOI: 10.1038/s41588-021-00911-1   
# Prepared for Statistical Genomics Workshop
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

# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)


##############################################################
# 1. Load the Visium dataset
# This folder must contain:
#   ├── filtered_feature_bc_matrix.h5
#   └── spatial/
##############################################################

data_dir <- "Module_4.3/Data/spatial/"

# Load the Visium dataset using the .h5 file + spatial image folder
breast <- Load10X_Spatial(
  data.dir = data_dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "breast"
)

##############################################################
# 2. Explore the raw data – Understanding what a Visium object contains
##############################################################

# Basic structure of the Seurat object
breast

# View the first few rows of the spatial metadata
head(breast@meta.data)

# Plot QC metrics: number of UMIs per spot
VlnPlot(breast, features = c("nCount_Spatial", "nFeature_Spatial"))

# Spatial visualization of nCount_Spatial on the tissue image
plot2 <- SpatialFeaturePlot(
  breast,
  features = "nCount_Spatial", 
  pt.size.factor = 3
) + theme(legend.position = "right")
plot2

##############################################################
# 3. SCTransform normalization
# This step normalizes the data while preserving spatial structure.
##############################################################

breast <- SCTransform(
  breast,
  assay = "Spatial",
  verbose = FALSE
)


##############################################################
# 4. Explore expression of selected marker genes
##############################################################

# Example gene (replace with breast-cancer relevant markers if needed)
SpatialFeaturePlot(breast, features = c("ERBB2", "MUC1", "KRT8"), pt.size.factor = 3)

##############################################################
# 5. High-quality feature plot
##############################################################

plot <- SpatialFeaturePlot(
  breast,
  features = "KRT8", pt.size.factor = 3
) + theme(
  legend.text  = element_text(size = 8),
  legend.title = element_text(size = 12),
  legend.key.size = unit(0.5, "cm")
)

jpeg(
  filename = "Module_4.3/Results/spatial_KRT8_example.jpg",
  height = 700,
  width = 1200,
  quality = 80
)
print(plot)
dev.off()

##############################################################
# 7. Dimensionality reduction + clustering
# Standard Seurat pipeline applied to spatial data
##############################################################

# PCA
breast <- RunPCA(breast, assay = "SCT", verbose = FALSE)

# Nearest neighbors + clustering
breast <- FindNeighbors(breast, reduction = "pca", dims = 1:30)
breast <- FindClusters(breast, verbose = FALSE, resolution = 2)

# UMAP
breast <- RunUMAP(breast, reduction = "pca", dims = 1:30)

##############################################################
# 8. Visualize clusters
##############################################################

p1 <- DimPlot(breast, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(breast, label = TRUE, label.size = 3, pt.size.factor = 2.7)

p1 + p2

# Highlight selected clusters on tissue
SpatialDimPlot(
  breast,
  cells.highlight = CellsByIdentities(breast, idents = c(2,1,4,3,5,8)),
  facet.highlight = TRUE,
  ncol = 3,
  pt.size.factor = 3
)

#check whether there are any specific clusters you are interested in!

##############################################################
# 9. Differential expression between clusters
##############################################################

de_markers <- FindMarkers(breast, ident.1 = 5, ident.2 = 6)
head(de_markers)

# Visualize top spatial markers
SpatialFeaturePlot(
  breast,
  features = rownames(de_markers)[1:3],
  alpha = c(0.1, 1),
  ncol = 3,
  pt.size.factor = 2.7
)

##############################################################
# 10. Identify spatially variable features
##############################################################

breast <- FindSpatiallyVariableFeatures(
  breast,
  assay = "SCT",
  features = VariableFeatures(breast)[1:1000],
  selection.method = "moransi"
)

# View top 6 spatially variable genes
top.features <- head(SpatiallyVariableFeatures(breast, selection.method = "moransi"), 6)
top.features

# Plot spatially variable genes
SpatialFeaturePlot(
  breast,
  features = top.features,
  ncol = 3,
  pt.size.factor = 3
)

