#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-LZ-1373
## Authors: Ayushi Agrawal
##
## Script Goal: Normalize using SCTransform, use input PCA dim and resolution
## to cluster cells and use ScType to assign cell types to the clusters.
##
## Usage example:
## Rscript 05_02a_subclustering_cell_type.R
##  --input 'input_seurat_object.RDS' \  # Input Seurat object 
##  --output '/output_directory' \       # Location for output files
##  --output_prefix "outputs_scRNA_" \   # Prefix for output files
##  --tissue "Brain" \                   # Tissue type to supply to ScType
##  --ndim 20 \                          # Number of PCs to use
##  --resolution 0.02 \                  # Resolution to use for Seurat clustering
##  --batch_var \                        # Batch variable for Harmony correction
##  --custom_db "custom_db.xlsx"         # Path to custom database for ScType
##
## Run "Rscript 05_02a_subclustering_cell_type.R --help" for more information
###############################################################################

# Get input arguments -----------------------------------------------------
library(optparse)
option_list <- list(
  make_option(c("-i", "--input"),
              action = "store", default = NA, type = "character",
              help = "Input Seurat object in RDS format (required)"
  ),
  make_option(c("-o", "--output"),
              action = "store", default = NA, type = "character",
              help = "Output directory, will create if it doesn't exist (required)"
  ),
  make_option(c("--output_prefix"),
              action = "store",
              default = "clustering_cell_type",
              type = "character",
              help = "Prefix for output files, [default %default]"
  ),
  make_option(c("--tissue"),
              action = "store", default = NA, type = "character",
              help = "Tissue type to supply to ScType (required)"
  ),
  make_option(c("-n", "--ndim"),
              action = "store", default = NA, type = "numeric",
              help = "Number of principal components to use (required)"
  ),
  make_option(c("-r", "--resolution"),
              action = "store", default = NA, type = "numeric",
              help = "Resolution parameter for FindClusters (required)"
  ),
  make_option("--batch_var",
              action = "store", default = NA, type = "character",
              help = "If provided, the batch variable will be used as the grouping variable for harmony batch correction, sample_label must also be provided (optional)"
  ),
  make_option("--custom_db",
              action = "store", default = NA, type = "character",
              help = "Path to custom database for ScType, must be a xlsx in the same format as the original (optional)"
  ),
  make_option("--remove_cluster",
              action = "store", default = NA, type = "numeric",
              help = "Seurat cluster number to remove"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Check if required args are provided
if (is.na(opt$resolution) | is.na(opt$input) | is.na(opt$output) | is.na(opt$ndim) | is.na(opt$tissue)) {
  stop("Missing one or more required arguments")
}

# Load required packages
library(dplyr)
library(Seurat)
library(HGNChelper)
library(openxlsx)
library(magrittr)
library(ggplot2)
library(readr)

# create the results folders
if (!(dir.exists(opt$output))) {
  dir.create(opt$output,recursive = TRUE)
}

#set seed
set.seed(42)

#increase RAM usage so the futurized Seurat functions... 
#can access larger global variables
options(future.globals.maxSize = 4000 * 1024^2)


# Cluster using optimized resolution --------------------------------------

# Read in the Seurat object
data <- readRDS(opt$input)

# Remove any previous normalization
DefaultAssay(data) <- "RNA"
assays_present <- Assays(data)
assays_keep <- setdiff(assays_present, "SCT")
data <- DietSeurat(data, assay = assays_keep)

if(!is.na(opt$remove_cluster)){
  print(table(data$seurat_clusters))
  print(paste0("Removing Seurat cluster ", opt$remove_cluster))
  data <- subset(data, subset = seurat_clusters==opt$remove_cluster, invert = TRUE)
  print(table(data$seurat_clusters))
}

if (!is.na(opt$batch_var)) {
  library(harmony)
  print(paste0("***** Correcting for batch effects using ", opt$batch_var, " *****"))
  data <- SCTransform(data, vst.flavor = "v2")
  data <- RunPCA(data,npcs = opt$ndim) %>%
    RunHarmony(assay.use="SCT",
               group.by.vars = opt$batch_var,
               kmeans_init_nstart=20,
               kmeans_init_iter_max=100) %>%
    RunUMAP(reduction = "harmony", dims = 1:opt$ndim)
  data <- FindNeighbors(
    object = data,
    dims = 1:opt$ndim,
    reduction = "harmony"
  )      
  data <- FindClusters(object = data, resolution = opt$resolution)
  print("***** Clustering and Normalization completed! (Harmony reduction) *****")
  
} else {
  print("***** No batch variable provided, skipping batch correction *****")
  
  
  data <- SCTransform(data, vst.flavor = "v2", verbose = FALSE)
  print("***** SCT normalization completed! *****")
  data <- RunPCA(data, npcs = opt$ndim, assay = "SCT")
  data <- FindNeighbors(
    object = data,
    dims = 1:opt$ndim,
    verbose = FALSE
  )
  data <- FindClusters(object = data, resolution = opt$resolution)
  data <- RunUMAP(object = data, dims = 1:opt$ndim)
  print("***** Clustering completed! *****")
}



# Run ScType --------------------------------------------------------------

source("/opt/sc-type/R/gene_sets_prepare.R")
source("/opt/sc-type/R/sctype_score_.R")

# If a custom database is provided, use that, otherwise use the default
if (!is.na(opt$custom_db)) {
  gs_list <- gene_sets_prepare(opt$custom_db, opt$tissue)
} else {
  gs_list <- gene_sets_prepare("/opt/sc-type/ScTypeDB_full.xlsx", opt$tissue)
}

# get cell-type by cell matrix
es.max <- sctype_score(
  scRNAseqData = as.matrix(data[["SCT"]]@counts), scaled = FALSE,
  gs = gs_list$gs_positive, gs2 = gs_list$gs_negative
)
rm(gs_list)
print("---ScType Score completed---")
# merge by cluster
cL_resutls <- do.call(
  "rbind",
  lapply(
    unique(data@meta.data$seurat_clusters),
    function(cl) {
      es.max.cl <- sort(rowSums(es.max[, rownames(data@meta.data[data@meta.data$seurat_clusters == cl, ])]), decreasing = !0)
      head(
        data.frame(
          cluster = cl,
          type = names(es.max.cl),
          scores = es.max.cl,
          ncells = sum(data@meta.data$seurat_clusters == cl)
        ),
        10
      )
    }
  )
)
# clear memory
rm(es.max)

sctype_scores <- cL_resutls %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = scores)

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells / 4] <- "Unknown"
print(sctype_scores[, 1:3])

write_csv(sctype_scores[order(sctype_scores$cluster),],
          file = paste0(
            opt$output,
            "/",
            opt$output_prefix,
            "_sctype_table.csv"
          )
)

data@meta.data$cell_type_sctype <- ""
for (j in unique(sctype_scores$cluster)) {
  cl_type <- sctype_scores[sctype_scores$cluster == j, ]
  data@meta.data$cell_type_sctype[data@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
}

# adding additional metadata of interest
data@meta.data$seurat_cluster_celltype_sctype <- paste0(data@meta.data$seurat_clusters, 
                                                        "_", data@meta.data$cell_type_sctype)
# Extract the numeric values from the strings
numeric_values <- as.numeric(gsub("[^0-9.]", "", unique(data@meta.data$seurat_cluster_celltype_sctype)))
# Sort the numeric values
sorted_vec <- unique(data@meta.data$seurat_cluster_celltype_sctype)[order(numeric_values)]
# Convert the metadata column to a factor with levels set accordingly
data@meta.data$seurat_cluster_celltype_sctype <- factor(data@meta.data$seurat_cluster_celltype_sctype, 
                                                        levels = sorted_vec)


rm(sctype_scores)
saveRDS(data,
        file = paste0(
          opt$output,
          "/",
          opt$output_prefix,
          "_clustered_and_cell_typed.rds"
        )
)




# cells per sample per cluster -------------------------------------------------
smp_cluster_counts <- unique(data@meta.data %>%
              group_by(sample_name) %>%
              mutate(total_numbers_of_cells_per_sample = n()) %>%
              group_by(seurat_clusters, .add=TRUE) %>%
              mutate(number_of_cells_per_sample_in_cluster = n(), 
                     proportion=
                       number_of_cells_per_sample_in_cluster/total_numbers_of_cells_per_sample) %>%
              select(sample_name,
                     library_prep_batch,
                     seurat_clusters,
                     total_numbers_of_cells_per_sample,
                     number_of_cells_per_sample_in_cluster, 
                     proportion))

pdf(file.path(opt$output, paste0(opt$output_prefix, "_proportion_of_cells_per_sample_cluster.pdf")),
    width = 25,
    height = 7)
print(ggplot(smp_cluster_counts, aes(x=seurat_clusters, y = round(proportion*100,1), fill = sample_name)) + 
        geom_bar(stat="identity", position = position_dodge(width = 0.8), width = 0.6, color = "black")+
        labs(x = "Seurat Cluster", y = "Cell Proportion", fill = "Sample") +
        theme_bw())
dev.off()




# Visualize the metadata columns that not numeric ------------------------------
meta_cols <- names(data@meta.data)[!sapply(data@meta.data, is.numeric)]
# Exclude "orig.ident"
meta_cols <- meta_cols[meta_cols != "orig.ident"]
# Exclude "DF.classifications_" columns
# meta_cols <- meta_cols[!grepl("^DF\\.classifications_", meta_cols)]
for (meta in meta_cols) {
  # Skip if there is only one value
  if (length(unique(data@meta.data[[meta]])) == 1){
    next
  }
  print(paste0("***** Plotting ", meta, " *****"))
  pdf(file.path(opt$output, paste0(opt$output_prefix, "_", meta, "_UMAP.pdf")),
      width = 12,
      height = 7)
  print(DimPlot(data,
                raster = FALSE,
                order = TRUE,
                label = FALSE,
                group.by = meta,
                reduction = "umap"
  ))
  dev.off()
  h <- ifelse(length(unique(data@meta.data[[meta]])) == 2,
              7,
              3 * ceiling(length(unique(data@meta.data[[meta]]))/2))
  pdf(file.path(opt$output, paste0(opt$output_prefix, "_", meta, "_split_UMAP.pdf")),
      width = 12,
      height = h)
  print(DimPlot(data,
                raster = FALSE,
                order = TRUE,
                label = FALSE,
                group.by = meta,
                split.by = meta,
                reduction = "umap",
                ncol = 2
  ) + NoLegend())
  dev.off()
  
}




# Visualize the metadata columns that are numeric ------------------------------
# Feature plots for numeric data
for(meta_feature in c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","nCount_SCT","nFeature_SCT",
                      grep("^SCT_norm_expr", colnames(data@meta.data), value = TRUE))){
  FeaturePlot(data,
              raster = FALSE,
              order = TRUE,
              label = FALSE,
              features = meta_feature,
              reduction = "umap") %>%
    ggsave(file = file.path(opt$output,
                            paste0(opt$output_prefix, "_", meta_feature, "_featureplot.pdf")),
           plot = .,
           width = 12,
           height = 7)
}




# UMAP with cell type labels --------------------------------------------------------------
pdf(file.path(opt$output, paste0(opt$output_prefix, "_seurat_clusters_labels_UMAP.pdf")),
    width = 12,
    height = 7
)
print(DimPlot(data,
              group.by = "seurat_clusters",
              raster = FALSE,
              order = TRUE,
              label = TRUE,
              reduction = "umap"
))
dev.off()

pdf(file.path(opt$output, paste0(opt$output_prefix, "_seurat_clusters_celltype_sctype_lables_UMAP.pdf")),
    width = 12,
    height = 7
)
print(DimPlot(data,
              group.by = "seurat_cluster_celltype_sctype",
              raster = FALSE,
              order = TRUE,
              label = TRUE,
              reduction = "umap"
))
dev.off()




# Find Markers --------------------------------------------------------------

print("**** FindAllMarkers start ****")
data <- PrepSCTFindMarkers(data)
# Find differentially expressed genes for identifying the subclusters
marker_res <- FindAllMarkers(data, 
                             assay = "SCT", 
                             slot = "data", 
                             test.use = "wilcox",
                             min.pct = 0.25)
marker_res %>%
  write.csv(., paste0(opt$output, "/", opt$output_prefix, "_FindAllMarkers_marker_genes_per_cluster.csv"), 
            row.names = FALSE)
print("**** FindAllMarkers end ****")





# Write session info
writeLines(
  capture.output(sessionInfo()),
  file.path(opt$output, "sessionInfo.txt")
)
print("***** Script Complete *****")


# END -----------------------------------------------------------