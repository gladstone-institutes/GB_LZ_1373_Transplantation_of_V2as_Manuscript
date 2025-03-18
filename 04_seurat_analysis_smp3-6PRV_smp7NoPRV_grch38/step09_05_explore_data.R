# this script was run locally on Ayushi's laptop
setwd("~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/")


library(dplyr)
library(Seurat)
library(HGNChelper)
library(magrittr)
library(ggplot2)
library(readr)
library(harmony)


# create the results folders
outdir <- "05_explore_data"
if (!(dir.exists(outdir))) {
  dir.create(outdir,recursive = TRUE)
}

# read in the object and marker list
markers <- read.csv("04_clustering_cell_type_no_cellbender_no_doubletfinder/20PC_0.6res/smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_FindAllMarkers_marker_genes_per_cluster.csv")
markers_rmBatch <- read.csv("04_clustering_cell_type_no_cellbender_no_doubletfinder/20PC_0.6res_rmBatchEffect/smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_rmBatchEffect_FindAllMarkers_marker_genes_per_cluster.csv")
dat <- readRDS("04_clustering_cell_type_no_cellbender_no_doubletfinder/20PC_0.6res/smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_clustered_and_cell_typed.rds")


# check if any mitochondiral genes? --------------------------------------------
# heatmap for top DE markers per cluster
markers %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
pdf(file.path(outdir, "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_20PC_0.6res_heatmap.pdf"),
    height = 20,
    width = 22)
print(DoHeatmap(dat, features = top10$gene) + NoLegend())
dev.off()


# Look for human cells only ----------------------------------------------------
# Subset by Y chromsome
human_cells = WhichCells(dat, slot = 'counts', expression = USP9Y > 0)

human_cells_UTY = WhichCells(dat, slot = 'counts', expression = UTY > 0)

human_cells_both = WhichCells(dat, slot = 'counts', expression = USP9Y > 0 & UTY > 0)

human_cells_either = WhichCells(dat, slot = 'counts', expression = USP9Y > 0 | UTY > 0)

data_human_USP9Y = subset(dat, cells = human_cells)
data_human_UTY = subset(dat, cells = human_cells_UTY)
data_human_both <- subset(dat, cells = human_cells_both)

pdf(file.path(outdir, "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_20PC_0.6res_Ychr_genes_violinplot_SCTdata.pdf"),
    width = 10)
print(VlnPlot(dat, features = c("USP9Y","UTY", "RPS4Y1"), group.by = "sample_name"))
dev.off()

pdf(file.path(outdir, "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_20PC_0.6res_Ychr_genes_violinplot_RNAcounts.pdf"),
    width = 10)
print(VlnPlot(dat, features = c("USP9Y","UTY", "RPS4Y1"), group.by = "sample_name",
              assay = "RNA", slot="counts"))
dev.off()


# Check for biological/batch effects  ------------------------------------------
dat_rmBatchEffect <- readRDS("04_clustering_cell_type_no_cellbender_no_doubletfinder/20PC_0.6res_rmBatchEffect/smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_rmBatchEffect_clustered_and_cell_typed.rds")

type_no <- dat@meta.data[,c("seurat_clusters","cell_type_sctype","seurat_cluster_celltype_sctype")]
type_no$cell_id <- rownames(type_no) 
type_rmBatchEffect <- dat_rmBatchEffect@meta.data[,c("seurat_clusters","cell_type_sctype","seurat_cluster_celltype_sctype")]
type_rmBatchEffect$cell_id <- rownames(type_rmBatchEffect) 
all <- merge(type_no, type_rmBatchEffect, by = "cell_id")

library(ggplot2)
library(reshape2)

# Transform the matrix in long format
data = all %>% group_by(seurat_clusters.x,seurat_clusters.y) %>% 
  summarise(n=log10(n())) %>% as.data.frame()

# plot the data
pdf(file.path(outdir, "cluster_membership_pre_vs_post_batch_correction.pdf"),
    width = 18,
    height = 14)
print(ggplot(data, aes(x = seurat_clusters.x, y = seurat_clusters.y, fill = n)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(10^n,2)), color = "white", size = 4) +
  scale_fill_gradient(low = "#33BBEE", high = "#EE3377") +
  ggtitle("Number of cells in the same cluster before and after batch correction") +
  xlab("Cluster labels before batch correction") + ylab("Cluster labels after batch correction") +
  theme_bw())
dev.off()


# UL transcript
pdf(file.path(outdir, "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_20PC_0.6res_ULgenes_violinplot.pdf"),
    width = 14)
print(VlnPlot(dat, features = c("UL18","UL19","UL25"), group.by = "sample_name"))
dev.off()



# correct for library
dat_lib_correct <- dat %>% RunHarmony(assay.use="SCT",
                                      group.by.vars = "library_prep_batch",
                                      kmeans_init_nstart=20,
                                      kmeans_init_iter_max=100) %>%
  RunUMAP(reduction = "harmony", dims = 1:20)

dat_lib_correct$lib_smp <- paste0(dat_lib_correct$library_prep_batch, dat_lib_correct$sample_name)
saveRDS(dat_lib_correct,
        file.path(outdir, "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_20PC_0.6res_rmBatchEffect_library.rds"))

pdf(file.path(outdir, "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_20PC_0.6res_rmBatchEffect_library_UMAP.pdf"),
    width = 14)
print(DimPlot(dat_lib_correct, group.by = "library_prep_batch"))
print(DimPlot(dat_lib_correct, group.by = "library_prep_batch", split.by = "library_prep_batch"))
print(DimPlot(dat_lib_correct, group.by = "sample_name", split.by = "sample_name"))
print(DimPlot(dat_lib_correct, group.by = "lib_smp"))
print(DimPlot(dat_lib_correct, split.by =  "sample_name"))
dev.off()


# Write session info
writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "sessionInfo.txt")
)
print("***** Script Complete *****")



# END -----------------------------------------------------------
