#!/usr/bin/env Rscript
#####################################################################################################
## Project ID: GB-LZ-1373
## Authors: Ayushi Agrawal
##
## Script Goal: Merge and normalize all Seurat objects; and visualize metadata for batch effects
##
## Usage example:
## Rscript 04_02_merge_and_visualize.R \
## --input_dir 01_qc/ \                                           # Directory containing the QC filtered Seurat objects
## --output_dir 02_merge_and_visualize \                          # Output directory
## --output_prefix "samples_3-6_grch38" \                         # Prefix for all output data and files
## --project "LZ_1373" \                                          # Project ID
## --exclude_samples LZ4087_07A_human.rds,LZ4087_07B_human.rds \  # Comma separated list of samples to exclude (optional)
## --npcs 15                                                      # Number of PCs to use for UMAP (optional)
## 
## Run "Rscript 04_02_merge_and_visualize.R --help" for more information
#####################################################################################################

# Get input arguments -------------------------------------------------------------------------------
library(optparse)
option_list <- list(
  make_option("--input_dir",
              action = "store", default = NA, type = "character",
              help = "Input directory (required)"
  ),
  make_option("--output_dir",
              action = "store", default = NA, type = "character",
              help = "Output directory (required)"
  ),
  make_option("--output_prefix",
              action = "store", default = NA, type = "character",
              help = "Prefix for all output data and files (required)"
  ),
  make_option("--project",
              action = "store", default = NA, type = "character",
              help = "Project ID for the Seurat objects (required)"
  ),
  make_option("--exclude_samples",
              action = "store", default = NA, type = "character",
              help = "Comma separated list of samples to exclude, the sample names should match the RDS file names (optional)"
  ),
  make_option("--npcs",
              action = "store", default = 20, type = "integer",
              help = "Number of PCs to use for UMAP, default is 20 (optional)"
  )
)
# Read in the arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Check that all required arguments are present
if (is.na(opt$input_dir) | is.na(opt$output_dir) | is.na(opt$output_prefix) | is.na(opt$project)) {
  stop("***** ERROR: Missing required arguments! *****")
}


# Load required packages
library(Seurat)
library(dplyr)
library(ggplot2)

# Set the working directory
setwd(opt$input_dir)

# Create the output directory if it doesn't exist
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}




# Merge the data into a single Seurat object --------------------------------------------------------------
# Read in all the sample data sets
samples_list <- list.files(".", pattern = "\\.rds$", full.names = TRUE, recursive=TRUE)

# If the exclude_samples argument is provided, remove these samples from the samples_list
if (!is.na(opt$exclude_samples)) {
  temp <- strsplit(as.character(opt$exclude_samples), ',')
  exclude_samples <- as.character(unlist(temp))
  print(exclude_samples)
  samples_list <- samples_list[!(basename(samples_list) %in% exclude_samples)]
  rm(temp)
}

sample_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(samples_list))

data_list <- lapply(samples_list, readRDS)
names(data_list) <- sample_names 

# Merge all the data sets
data <- merge(data_list[[1]], 
              y = unlist(data_list[2:length(data_list)]), 
              add.cell.ids = sample_names,
              project = opt$project, 
              merge.data = FALSE)

rm(samples_list, sample_names, data_list)

# Save the processed dataset
saveRDS(data, 
        file = file.path(opt$output_dir, 
                         paste0(opt$output_prefix, 
                                "_merged_data.rds")))

print("***** Dataset merge completed! *****")




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
                                "_merged_data_sct_pca_umap.rds")))
print("***** SCT normalization, PCA and UMAP completed! *****")




# Visualize the data ------------------------------------------------------------------------------------
# PCA plots
pdf(file = file.path(opt$output_dir, paste0(opt$output_prefix, "_PCA_plots.pdf")),
    width = 12, height = 7)
print(ElbowPlot(data, ndims = 50))
print(ElbowPlot(data, ndims = 20))
print(VizDimLoadings(data, dims = 1:2, reduction = "pca"))
print(DimPlot(data, reduction = "pca", raster = FALSE))
print(DimHeatmap(data, dims = 1:2, cells = 500, balanced = TRUE))
print(DimHeatmap(data, dims = 3:4, cells = 500, balanced = TRUE))
dev.off()


# Quantitative approach to identify PCA plot elbow
# Reference: https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# Determine percent of variation associated with each PC
pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation 
# associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# Minimum of the two calculation
pcs <- min(co1, co2)
print(paste0("***** ", pcs+1, " PCs are recommended for this data using quantitative approach."))


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
                      "SCT_norm_expr_ChR2.H134R.", "SCT_norm_expr_mRFP1", "SCT_norm_expr_UL25", 
                      "SCT_norm_expr_UL21", "SCT_norm_expr_UL19", "SCT_norm_expr_UL18", 
                      "SCT_norm_expr_UL15", "SCT_norm_expr_YFP", "SCT_norm_expr_UL20")){
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




# save the session info
writeLines(capture.output(sessionInfo()), 
           file.path(opt$output_dir, "sessionInfo.txt"))

# record logs
print("***** Visualization completed! *****")


########################## END ##########################  
