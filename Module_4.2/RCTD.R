# Install R packages (only need to do this once, or if you need to update packages)
#install.packages("devtools")
#devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
#install.packages("Seurat")
#install.packages("ggplot2")
#install.packages("SeuratObject", type = "source")
#install.packages("hdf5r")

# Load libraries
library(spacexr)
library(Seurat)
library(SeuratObject)
library(ggplot2) #for saving result plots
library(Matrix)
library(dplyr)
library(hdf5r)

# Defining the spatial directory (replace with your paths)
visium_dir <- "/cloud/project/Module_4.3/Data/spatial"

# Reading the raw coordinates
coords_path <- file.path(visium_dir, "spatial", "tissue_positions_list.csv")

vis_coords <- read.csv(coords_path, header = FALSE)
colnames(vis_coords) <- c("barcode", "in_tissue", "array_row", "array_col",
                          "pxl_row_in_fullres", "pxl_col_in_fullres")

head(vis_coords)

#We need to flip teh coordinates to match the orientation of the tissue on the slide with spacexr result plots (it is optional)

# Multiply x and y coordinates by -1 to flip/mirror
vis_coords[,3] <- vis_coords[,3]*-1
vis_coords[,4] <- vis_coords[,4]*-1
head(vis_coords)
# the x-y coordinates are now flipped

# Resave file.
# Note that due to the file name requirement, we will write over the original file with this new one
write.table(vis_coords,
            coords_path,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = ",")

# Load Visium data 
VisiumData <- read.VisiumSpatialRNA("/cloud/project/Module_4.3/Data/spatial/")

# Extract barcodes
barcodes <- colnames(VisiumData@counts)

# Plot UMIs per spot for quality check 

plot_puck_continuous(
  puck = VisiumData,
  barcodes = barcodes,
  plot_val = VisiumData@nUMI,
  size = 1,
  ylimit = c(0, round(quantile(VisiumData@nUMI, 0.9))),
  title = "nUMI per spot"
)

# Now, we need to read the reference single cell RNA-seq dataset
Counts <- Read10X_h5("/cloud/project/Module_4.3/Data/scRNA/filtered_feature_bc_matrix.h5")

sc <- CreateSeuratObject(Counts)

head(sc@assays$RNA@layers$counts)
sc_counts <- sc@assays$RNA@layers$counts
# Seurat v5 stores gene names in @features
gene_names <- rownames(sc@assays$RNA@features)

# Assign to your count matrix
rownames(sc_counts) <- gene_names
# Assign barcodes from the Seurat object to the counts matrix
colnames(sc_counts) <- colnames(sc)
head(colnames(sc_counts))
# Now loading the metadata for the cell type annotation

meta <- read.csv("/cloud/project/Module_4.3/Data/scRNA/metadata.csv", row.names=1)

# Ensure barcodes match
meta <- meta[colnames(sc), , drop=FALSE]

# Add metadata
sc <- AddMetaData(sc, metadata = meta)

head(sc$celltype_major)

cell_types <- sc$celltype_major
head(cell_types)

# Check how many cells are present for each cell type. For RCTD, every cell type must have >= 25 cells
table(cell_types)

# Removing unsafe characters such as /, spaces, hyphens

cell_types_clean <- gsub("/", "_", cell_types)
cell_types_clean <- gsub("-", "_", cell_types_clean)
cell_types_clean <- gsub(" ", "_", cell_types_clean)

# Convert to factor
cell_types_clean <- as.factor(cell_types_clean)

# Preparing nUNI per barcode

sc_umis <- setNames(sc$nCount_RNA, colnames(sc))
head(sc_umis)

# Creating reference object for RCTD

SCreference <- Reference(
  sc_counts,        # dgCMatrix, genes x cells
  cell_types_clean, # factor with cell type per barcode
  sc_umis           # named vector of UMI counts
)

SCreference

# The processes below require more RAM and take quite some time to run
# So I have commented that out. You can try running it when you have time
# We will load the processed object directly 

# Create and run RCTD
#myRCTD <- create.RCTD(VisiumData, 
#                      SCreference,
#                      max_cores = 10      # adjust as needed
#                      )

# Run RCTD
#myRCTD <- run.RCTD(
#  myRCTD,
#  doublet_mode = "full"   # or "doublet" or "cell"
#)

# Loading the processed object 

myRCTD <- readRDS("Module_4.3/Data/RCTD_processed.rds")

# Creating output folder
resultsdir <- "/cloud/project/Module_4.3/Results/"

# Extracting the spatial barcodes and weights
barcodes <- colnames(myRCTD@spatialRNA@counts)
weights <- myRCTD@results$weights   # cell type proportions

# Normalizing the weights so each spot sums to 1
norm_weights <- normalize_weights(weights)

# Cell type names object
cell_type_names <- colnames(norm_weights)

# Inspecting two example spots (optional)
subset_df <- as.data.frame(t(as.data.frame(norm_weights[1:2,])))
subset_df$celltypes <- rownames(subset_df)
subset_df

# Creating plots of the cell type proportions plotted on the spatial spots and saving them in the results folder

for (i in seq_along(cell_type_names)) {
  
  ct <- cell_type_names[i]
  
  plot_puck_continuous(
    puck = myRCTD@spatialRNA,
    barcodes = barcodes,
    plot_val = norm_weights[, ct],
    title = ct,
    size = 1
  )
  
  ggsave(
    filename = paste0(resultsdir, ct, "_weights.jpg"),
    height = 5,
    width = 5,
    units = "in",
    dpi = 300
  )
}

# For each spot, find the cell type with the highest weight
dominant_ct <- apply(norm_weights, 1, function(x) {
  names(x)[which.max(x)]
})

# Add dominant cell type to RCTD object (optional)
myRCTD@results$dominant_type <- dominant_ct

#Storing teh coordinates for plotting
coords <- myRCTD@spatialRNA@coords
colnames(coords) <- c("x", "y")

# Dominant cell dataframe for plotting
df_dom <- data.frame(
  barcode = rownames(coords),
  x = coords$x,
  y = coords$y,
  dominant = dominant_ct
)


ggplot(df_dom, aes(x = x, y = y, color = dominant)) +
  geom_point(size = 1.5) +
  theme_void() +
  ggtitle("Dominant Cell Type per Visium Spot") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )


