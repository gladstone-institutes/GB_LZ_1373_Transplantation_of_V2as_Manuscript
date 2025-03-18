#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-LZ-1373
## Authors: Ayushi Agrawal
##
## Script Goal: Subset the snRNA-seq RDS object for select select clusters and
##          visualize in UMAP space to check for batch effects.
##
## Usage example:
## Rscript 10_01_subset_visualize.R
##  --input_processed 'seurat_object.RDS' \  # Input Seurat object 
##  --output '/output_directory' \           # Location for output files
##  --output_prefix "outputs_scRNA_" \       # Prefix for output files
##  --npcs 20 \                              # Number of PCs to use
##  --cluster_subset "2,3"                   # Cluster(s) to subset and run clustering for
##  --gene_subset "RBFOX3"                   # Gene >0 exp to subset and run clustering for
##
## Run "Rscript 10_01_subset_visualize.R --help" for more information
###############################################################################

# Get input arguments -----------------------------------------------------
library(optparse)
option_list <- list(
  make_option("--input_processed",
              action = "store", default = NA, type = "character",
              help = "Input Seurat object post clustering in RDS format (required)"
  ),
  make_option("--input_orig",
              action = "store", default = NA, type = "character",
              help = "Merged Seurat object with no processing in RDS format (required)"
  ),
  make_option("--output_dir",
              action = "store", default = NA, type = "character",
              help = "Output directory (required)"
  ),
  make_option("--output_prefix",
              action = "store", default = "seurat_subset", type = "character",
              help = "Prefix for all output data and files (optional)"
  ),
  make_option(c("-n", "--npcs"),
              action = "store", default = NA, type = "numeric",
              help = "Number of principal components to use (required)"
  ),
  make_option(c("--cluster_subset"),
              action = "store", default = NA, type = "character",
              help = "The input will be subset by these seurat clusters (required)"
  ),
  make_option(c("--gene_subset"),
              action = "store", default = NA, type = "character",
              help = "The input will be subset by the cells with >0 exp of this gene. (optional)"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Check if required args are provided
if (is.na(opt$input_processed) | is.na(opt$input_orig) | is.na(opt$output_dir) | is.na(opt$npcs) | is.na(opt$cluster_subset)) {
  stop("Missing one or more required arguments")
}




# Load required packages --------------------------------------------------
library(dplyr)
library(Seurat)
library(HGNChelper)
library(openxlsx)
library(magrittr)
library(ggplot2)
library(readr)

# create the results folders
if (!(dir.exists(opt$output_dir))) {
  dir.create(opt$output_dir,recursive = TRUE)
}

#set seed
set.seed(42)

#increase RAM usage so the futurized Seurat functions... 
#can access larger global variables
options(future.globals.maxSize = 4000 * 1024^2)




# Sub clustering using optimized resolution --------------------------------------

# Read in the Seurat object
data <- readRDS(opt$input_processed)
Idents(data) <- data$seurat_clusters
print("Seurat clusters for all data: ")
print(table(data$seurat_clusters))

# Read in unprocessed merged data
data_orig <- readRDS(opt$input_orig)

# Subset the Seurat object to include only the seurat clusters provided
temp <- strsplit(as.character(opt$cluster_subset), ',')
cell_type_clusters <- as.numeric(as.character(unlist(temp)))
cells_to_subcluster <- WhichCells(data, idents = cell_type_clusters)
data <- subset(x = data_orig, cells = cells_to_subcluster)
rm(data_orig, temp, cell_type_clusters, cells_to_subcluster)

if(!is.na(opt$gene_subset)){
  expr <- FetchData(object = data, vars = opt$gene_subset)
  data <- data[, which(x = rowSums(expr) > 0)]
}

print("***** Subsetting completed! *****")




# SCTransform normalization with UL transcripts --------------------------------------------------------------
data <- SCTransform(data, 
                    vst.flavor = "v2", 
                    verbose = TRUE, 
                    new.assay.name = "sctwithultranscripts") 

# Extract the expression of the UL transcripts
DefaultAssay(data) <- "sctwithultranscripts"
remove_genes <- c("ChR2(H134R)", "mRFP1", "UL15", "UL18", "UL19", "UL20", "UL21", "UL25", "YFP")

# Ensure the genes exist in the dataset before extracting
remove_genes <- intersect(rownames(data), remove_genes)
print("Genes to remove:")
print(remove_genes)

if (length(remove_genes) > 0) {
  expression_matrix <- GetAssayData(data, 
                                    assay = "sctwithultranscripts", 
                                    slot = "data")[remove_genes, , drop = FALSE]
  print("Expression matrix")
  print(dim(expression_matrix))
  print(head(expression_matrix))
  
  # Convert sparse matrix to a dense matrix
  expression_matrix <- as.matrix(expression_matrix)
  
  # Convert to a data frame and transpose so cells are rows
  expression_metadata <- as.data.frame(t(expression_matrix))
  
  # Rename columns to indicate they are expression values
  colnames(expression_metadata) <- paste0("SCT_norm_expr_", colnames(expression_metadata))
  
  # Add the expression values to metadata
  data <- AddMetaData(data, metadata = expression_metadata)
  
  print(head(data@meta.data))
  
} else {
  print("None of the remove_genes are present in the dataset.")
}




# Normalization, PCA, and UMAP after removing UL transcripts -------------------------------------------
# Remove the UL transcripts from the data 
DefaultAssay(data) <- "RNA"
data <- subset(data, 
               features = setdiff(rownames(data), remove_genes))

# Run SCTransform normalization, PCA and UMAP
data <- SCTransform(data, vst.flavor = "v2", verbose = TRUE, new.assay.name = "SCT") %>%
  RunPCA(assay = "SCT", npcs = 50) %>%
  RunUMAP(reduction = "pca", dims = 1:opt$npcs, verbose = FALSE)

print(data)

# Save the processed dataset
saveRDS(data, 
        file = file.path(opt$output_dir, 
                         paste0(opt$output_prefix, 
                                "_data_sct_pca_umap.rds")))
print("***** SCT normalization, PCA and UMAP completed! *****")




# Visualize the data ------------------------------------------------------------------------------------
# PCA plot
pdf(file = file.path(opt$output_dir, paste0(opt$output_prefix, "_PCA_plots.pdf")),
    width = 12, height = 7)
print(ElbowPlot(data, ndims = 50))
print(ElbowPlot(data, ndims = 20))
print(VizDimLoadings(data, dims = 1:2, reduction = "pca"))
print(DimPlot(data, reduction = "pca", raster = FALSE))
print(DimHeatmap(data, dims = 1:2, cells = 500, balanced = TRUE))
print(DimHeatmap(data, dims = 3:4, cells = 500, balanced = TRUE))
dev.off()

# # Quantitative approach to identify PCA plot elbow
# # Reference: https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# # Determine percent of variation associated with each PC
# pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100
# # Calculate cumulative percents for each PC
# cumu <- cumsum(pct)
# # Determine which PC exhibits cumulative percent greater than 90% and % variation 
# # associated with the PC as less than 5
# co1 <- which(cumu > 90 & pct < 5)[1]
# # Determine the difference between variation of PC and subsequent PC
# co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# # Minimum of the two calculation
# pcs <- min(co1, co2)
# print(paste0(pcs+1, " PCs are recommended for this data using quantitative approach."))



# Function to generate visualizations for the metadata variables
plot_metadata <- function(sc_dat, metadata_variable){
  # Skip if there is only one value
  if (length(unique(sc_dat@meta.data[[metadata_variable]])) > 1){
    print(paste0("***** Plotting ", metadata_variable, " *****"))
    DimPlot(sc_dat, 
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            group.by=metadata_variable, 
            reduction = "umap") %>%
      ggsave(file = file.path(opt$output_dir, 
                              paste0(opt$output_prefix, "_", metadata_variable, "_umap.pdf")),
             plot = .,
             width = 12,
             height = 7)
    pdf_width <- (min(length(unique(sc_dat@meta.data[[metadata_variable]])),4)) * 6  
    pdf_height <- ifelse(length(unique(sc_dat@meta.data[[metadata_variable]])) <= 4,
                         7,
                         5 * ceiling(length(unique(sc_dat@meta.data[[metadata_variable]]))/4))
    (DimPlot(sc_dat, 
             raster = FALSE, 
             order = TRUE, 
             label = FALSE, 
             group.by=metadata_variable,
             split.by=metadata_variable, 
             reduction = "umap",
             ncol = 4) + NoLegend()) %>%
      ggsave(file = file.path(opt$output_dir, 
                              paste0(opt$output_prefix, "_", metadata_variable, "_split_umap.pdf")),
             plot = .,
             width = pdf_width,
             height = pdf_height,
             limitsize = FALSE)
  }
}

# UMAPs for the metadata variables
# get the metadata columns that not numeric
meta_cols <- names(data@meta.data)[!sapply(data@meta.data, is.numeric)]
# Exclude "orig.ident"
meta_cols <- meta_cols[meta_cols != "orig.ident"]
data@meta.data[meta_cols] <- lapply(data@meta.data[meta_cols], as.factor)
lapply(meta_cols, function(col_name) {
  plot_metadata(data, col_name)
})


# Feature plots for numeric data
for(meta_feature in c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","nCount_SCT","nFeature_SCT",
                      grep("^SCT_norm_expr", colnames(data@meta.data), value = TRUE))){
  FeaturePlot(data, 
              raster = FALSE, 
              order = TRUE, 
              label = FALSE, 
              features = meta_feature, 
              reduction = "umap") %>%
    ggsave(file = file.path(opt$output_dir, 
                            paste0(opt$output_prefix, "_", meta_feature, "_featureplot.pdf")),
           plot = .,
           width = 12,
           height = 7)
}



# Write session info
writeLines(
  capture.output(sessionInfo()),
  file.path(opt$output_dir, "sessionInfo.txt")
)
print("***** Script Complete *****")



# END -----------------------------------------------------------
