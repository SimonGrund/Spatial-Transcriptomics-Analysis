# Spatial Transcriptomics Analysis

This repository contains R scripts and data for analyzing single-cell RNA sequencing (scRNA-seq) and spatial transcriptomics data, developed for the Statistical Genomics Workshop. The project consists of two separate but related modules focused on breast cancer gene expression analysis.

## Overview

- **Module 4.1**: Single-Cell RNA-Seq Analysis using Seurat
- **Module 4.2**: Spatial Transcriptomics Analysis with Visium data

Both modules utilize breast cancer datasets and demonstrate complementary approaches to understanding tumor heterogeneity and spatial organization.

---

## Module 4.1: Single-Cell RNA-Seq Analysis

### Description

This module provides a comprehensive tutorial for analyzing single-cell RNA-sequencing data using the Seurat package. The workflow guides users through standard scRNA-seq analysis steps from quality control to cell type identification and marker gene discovery.

### Dataset

- **Input Files**:
  - `filtered_feature_bc_matrix.h5`: 10x Genomics format gene expression matrix
  - `metadata.csv`: Cell metadata including cell type annotations

### Key Analysis Steps

1. **Data Loading and Inspection**
   - Load 10x Genomics .h5 matrix format
   - Import and integrate cell metadata
   - Assess data sparsity

2. **Quality Control**
   - Calculate mitochondrial gene percentages
   - Visualize QC metrics (nFeature_RNA, nCount_RNA, percent.mt)
   - Optional filtering of low-quality cells

3. **Normalization and Feature Selection**
   - LogNormalize data with scale factor 10,000
   - Identify 2,000 highly variable features

4. **Dimensionality Reduction**
   - Scale data and run PCA
   - Generate elbow plot for PC selection
   - Visualize data in PCA space

5. **Clustering and UMAP**
   - Construct nearest neighbor graph (15 PCs)
   - Cluster at multiple resolutions (0.5 and 1.0)
   - Generate UMAP embeddings

6. **Cell Type Analysis**
   - Analyze cell type frequencies
   - Filter rare cell types (< 100 cells)
   - Create filtered cell type annotations

7. **Marker Gene Identification**
   - Identify markers for all cell types
   - Extract top markers per cluster
   - Visualize markers with DotPlot

### Dependencies

```r
library(Seurat)
```

### Usage

```r
# Navigate to project directory and run:
source("Module 4.1/scRNA_tutorial_code.R")
```

### Output

- QC plots (violin plots)
- PCA and elbow plots
- UMAP visualizations with cell type annotations
- Marker gene dot plots
- Cell type frequency tables

---

## Module 4.2: Spatial Transcriptomics Analysis

### Description

This module demonstrates spatial transcriptomics analysis using 10x Genomics Visium technology on a breast cancer sample (CID44971). The analysis includes both Seurat-based spatial analysis and cell type deconvolution using RCTD (Robust Cell Type Decomposition).

### Dataset

**Source**: Wu et al., Nature Genetics 2021 (DOI: [10.1038/s41588-021-00911-1](https://doi.org/10.1038/s41588-021-00911-1))

- **Spatial Data** (`Data/spatial/`):
  - `filtered_feature_bc_matrix.h5`: Gene expression counts per spot
  - `spatial/tissue_positions_list.csv`: Spatial coordinates
  - `spatial/tissue_hires_image.png`: High-resolution H&E image
  - `spatial/tissue_lowres_image.png`: Low-resolution H&E image
  - `spatial/scalefactors_json.json`: Image scale factors

- **Single-Cell Reference** (`Data/scRNA/`):
  - `filtered_feature_bc_matrix.h5`: scRNA-seq reference data
  - `metadata.csv`: Cell type annotations for reference

### Scripts

1. **`Seurat_spatial_basic.R`**
   - Basic spatial data loading and visualization
   - SCTransform normalization
   - Marker gene expression plotting

2. **`Spatial_Seurat_Final.R`**
   - Complete spatial analysis workflow with detailed QC
   - Dimensionality reduction and clustering
   - Differential expression analysis
   - Spatially variable gene detection using Moran's I
   - Interactive questions for workshop participants

3. **`RCTD.R`**
   - Cell type deconvolution using spacexr package
   - Integration of scRNA-seq reference with spatial data
   - Generation of cell type proportion maps
   - Dominant cell type assignment per spot

### Key Analysis Steps

#### Seurat Spatial Analysis

1. **Data Loading**
   - Load Visium dataset with spatial coordinates and images
   - Quality control metric visualization

2. **Quality Control**
   - Calculate mitochondrial percentage and capture efficiency
   - Spatial QC visualization on tissue
   - Optional filtering based on quantile thresholds

3. **Normalization**
   - SCTransform normalization (preserves spatial structure)

4. **Marker Gene Expression**
   - Visualize breast cancer markers (ERBB2, MUC1, KRT8)
   - High-quality spatial feature plots

5. **Clustering**
   - PCA-based dimensionality reduction
   - Graph-based clustering (resolution 2.0)
   - UMAP generation for visualization

6. **Differential Expression**
   - Compare gene expression between clusters
   - Identify cluster-specific markers

7. **Spatial Analysis**
   - Find spatially variable features using Moran's I
   - Identify genes with spatial autocorrelation

#### RCTD Deconvolution

1. **Data Preparation**
   - Coordinate transformation and orientation correction
   - Load Visium spatial data
   - Prepare scRNA-seq reference (cell types, counts, UMI counts)

2. **RCTD Execution**
   - Create RCTD object with spatial and reference data
   - Run deconvolution in "full" doublet mode
   - Normalize cell type weights per spot

3. **Visualization**
   - Generate cell type proportion maps for each cell type:
     - B_cells
     - CAFs (Cancer-Associated Fibroblasts)
     - Cancer_Epithelial
     - Endothelial
     - Myeloid
     - Normal_Epithelial
     - Plasmablasts
     - PVL (Perivascular-Like)
     - T_cells
   - Identify dominant cell type per spot

### Dependencies

```r
# Seurat analysis
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ape)
library(pheatmap)

# RCTD analysis
library(spacexr)
library(SeuratObject)
library(Matrix)
library(hdf5r)
```

### Installation

```r
# Install required packages (run once)
install.packages("devtools")
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
install.packages("Seurat")
install.packages("ggplot2")
install.packages("SeuratObject", type = "source")
install.packages("hdf5r")
install.packages("ape")
install.packages("pheatmap")
install.packages('BiocManager')
BiocManager::install('glmGamPoi')
```

### Usage

**Important**: Before running the scripts, you need to update the file paths to match your local directory structure:

1. **For Seurat_spatial_basic.R**: Change line 31 from `Module_4.3/Data/spatial/` to `Module_4.2/Data/spatial/` and line 95 from `Module_4.3/Results/` to `Module_4.2/Results/`

2. **For Spatial_Seurat_Final.R**: Uses correct `Module_4.2` paths and can be run as-is

3. **For RCTD.R**: Change all instances of `Module_4.3` to `Module_4.2` (lines 19, 48, 65, 81, 139, 142) and update absolute paths `/cloud/project/` to your working directory or use relative paths

Then run:

```r
# For Seurat spatial analysis (basic):
source("Module_4.2/Seurat_spatial_basic.R")

# For complete Seurat spatial analysis with questions:
source("Module_4.2/Spatial_Seurat_Final.R")

# For RCTD cell type deconvolution:
source("Module_4.2/RCTD.R")
```

**Note**: The RCTD script includes pre-processed results (`Data/RCTD_processed.rds`) to save computation time, as full RCTD analysis requires significant RAM and processing time.

### Output

- **Results Directory** (`Module_4.2/Results/`):
  - `spatial_KRT8_example.jpg`: KRT8 expression spatial plot
  - `*_weights.jpg`: Cell type proportion maps (10 cell types)
  - Quality control plots
  - Cluster visualizations
  - Spatially variable gene plots

- **HTML Report**:
  - `spatial_breast_cancer.html`: Compiled analysis report

---

## Relationship Between Modules

While both modules are separate analyses, they are complementary:

- **Module 4.1** provides the foundational skills for scRNA-seq analysis and cell type identification
- **Module 4.2** extends these concepts to spatial data, showing how cell types are organized in tissue space
- The scRNA-seq data in Module 4.2 serves as a reference for RCTD deconvolution, linking single-cell resolution to spatial spots

---

## Citation

If you use the spatial transcriptomics dataset, please cite:

Wu, S.Z., Al-Eryani, G., Roden, D.L. et al. A single-cell and spatially resolved atlas of human breast cancers. *Nat Genet* 53, 1334–1347 (2021). https://doi.org/10.1038/s41588-021-00911-1

---

## Requirements

- **R version**: 4.0 or higher recommended
- **RAM**: Minimum 8GB (16GB+ recommended for RCTD)
- **Disk space**: ~500MB for data and results

---

## Project Structure

```
.
├── Module 4.1/                             # Note: Directory name has a space
│   ├── filtered_feature_bc_matrix.h5
│   ├── metadata.csv
│   └── scRNA_tutorial_code.R
├── Module_4.2/                             # Note: Directory name has an underscore
│   ├── Data/
│   │   ├── spatial/
│   │   │   ├── filtered_feature_bc_matrix.h5
│   │   │   └── spatial/
│   │   │       ├── tissue_positions_list.csv
│   │   │       ├── tissue_hires_image.png
│   │   │       ├── tissue_lowres_image.png
│   │   │       ├── scalefactors_json.json
│   │   │       ├── detected_tissue_image.jpg
│   │   │       └── aligned_fiducials.jpg
│   │   ├── scRNA/
│   │   │   ├── filtered_feature_bc_matrix.h5
│   │   │   └── metadata.csv
│   │   └── RCTD_processed.rds
│   ├── Results/
│   │   ├── spatial_KRT8_example.jpg
│   │   └── [10 cell type weight maps]
│   ├── Seurat_spatial_basic.R
│   ├── Spatial_Seurat_Final.R
│   ├── RCTD.R
│   └── spatial_breast_cancer.html
└── README.md
```

---

## Contact

For questions or issues, please open an issue in the GitHub repository.

---

## License

This project is intended for educational purposes as part of a Statistical Genomics Workshop.
