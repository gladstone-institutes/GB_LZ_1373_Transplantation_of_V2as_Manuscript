#!/usr/bin/env Rscript

################################################################################
## Project ID: GB-LZ-1373
## Authors: Ayushi Agrawal
## Script Goal: Identify marker genes for neuronal mRFP+ve cells in subcluster 6
## Usage:Rscript 05_02b_split_neuronal_subcluster6.R
################################################################################


# Setup ------------------------------------------------------------------------
# load required libraries
library(Seurat)
library(HGNChelper)
library(dplyr)
library(openxlsx)
library(magrittr)
library(ggplot2)
library(readr)

# folder paths
basedir <- file.path("/gladstone/bioinformatics/projects",
                     "lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023",
                     "results/05_snRNAseq_subcluster_neurons_smp3-6PRV_smp7NoPRV")
setwd(file.path(basedir,
                "02_subclustering_cell_type/30PC_0.30res_rmBatchEffect"))

# create output directory
outdir <- "split_subcluster6"
if (!(dir.exists(outdir))) {
  dir.create(outdir,recursive = TRUE)
}

output_prefix <- "neurons_smp3-7_grch38_noCB_noDF_NoPRV_30PC_0.30res_rmBatchEffect_split_subcluster6"


# Main -------------------------------------------------------------------------
# read in the data 
data <- readRDS("neurons_smp3-7_grch38_noCB_noDF_NoPRV_30PC_0.30res_rmBatchEffect_clustered_and_cell_typed.rds")

# Create new metadata variable to track clusters
data$modified_clusters <- as.character(data$seurat_clusters)

# Identify cluster 6 cells
cluster6_cells <- which(data$seurat_clusters == "6")

# Subset cluster 6 cells and check the condition
data$modified_clusters[cluster6_cells] <- ifelse(data$SCT_norm_expr_mRFP1[cluster6_cells] > 0, "6_2", "6_1")

# Convert back to a factor if needed
cluster_levels <- c(as.character(0:5), "6_1", "6_2", as.character(7:(max(as.numeric(data$seurat_clusters)) - 1)))
data$modified_clusters <- factor(data$modified_clusters, levels = cluster_levels)

# Check the result
table(data$modified_clusters)

# Change the idents of the Seurat object
Idents(data) <- data$modified_clusters


# Run ScType --------------------------------------------------------------
source("/opt/sc-type/R/gene_sets_prepare.R")
source("/opt/sc-type/R/sctype_score_.R")

# If a custom database is provided, use that, otherwise use the default
gs_list <- gene_sets_prepare("/opt/sc-type/ScTypeDB_full.xlsx", "Brain")

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
            outdir,
            "/",
            output_prefix,
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
          outdir,
          "/",
          output_prefix,
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

pdf(file.path(outdir, paste0(output_prefix, "_proportion_of_cells_per_sample_cluster.pdf")),
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
  pdf(file.path(outdir, paste0(output_prefix, "_", meta, "_UMAP.pdf")),
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
  pdf(file.path(outdir, paste0(output_prefix, "_", meta, "_split_UMAP.pdf")),
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
    ggsave(file = file.path(outdir,
                            paste0(output_prefix, "_", meta_feature, "_featureplot.pdf")),
           plot = .,
           width = 12,
           height = 7)
}




# UMAP with cell type labels --------------------------------------------------------------
pdf(file.path(outdir, paste0(output_prefix, "_seurat_clusters_labels_UMAP.pdf")),
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

pdf(file.path(outdir, paste0(output_prefix, "_seurat_clusters_celltype_sctype_lables_UMAP.pdf")),
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
  write.csv(., paste0(outdir, "/", output_prefix, "_FindAllMarkers_marker_genes_per_cluster.csv"), 
            row.names = FALSE)

# heatmap
#heatmap for top 10 marker genes for each cluster
marker_res %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

pdf(file = paste0(outdir, "/", output_prefix, "_FindAllMarkers_marker_genes_per_cluster_heatmap.pdf"),
    width = 12,
    height = 15)
print(DoHeatmap(data, features = top10$gene) +
        scale_fill_gradientn(colors = c("#33BBEE", "white", "#EE3377")))
dev.off()

print("**** FindAllMarkers end ****")



# Write session info
writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "sessionInfo.txt")
)
print("***** Script Complete *****")


## END 