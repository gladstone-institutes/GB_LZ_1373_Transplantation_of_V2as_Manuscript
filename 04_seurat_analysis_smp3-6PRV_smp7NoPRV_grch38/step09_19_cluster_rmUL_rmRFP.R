#!/usr/bin/env Rscript
#run locally on Ayushi's laptop

#load required packages
library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)

#set seed
set.seed(42)

#increase RAM usage so the futurized Seurat functions... 
#can access larger global variables
options(future.globals.maxSize = 4000 * 1024^2)

# set working directory
setwd(paste0("~/Dropbox (Gladstone)/GB-LZ-1373/results/",
             "09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/",
             "09_clustering_cell_type_no_CB_no_DF_HumanCells/",
             "20PC_0.6res_rmBatchEffect_rmCluster16/"))

# 1. Take out the UL (UL18, 19, 25 and mRFP) transcripts from the total Neurons+glia 
#    snRNAseq dataset and see if this “Cluster 15” remains a unique cluster, or if the 
#    cells are now dispersed throughout the dataset -----------------------------------------------------
# visualize cluster 15
data <- readRDS(paste0("smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_",
                       "rmBatchEffect_rmCluster16_clustered_and_cell_typed.rds"))
data$clus15 <- ifelse(data$seurat_clusters == 15, "Yes", "No")
pdf("Cluster_15_in_full_data_withULproteins.pdf",
    width = 12,
    height = 7)
print(DimPlot(data, 
              group.by = "clus15", 
              raster = F,
              order = T,
              cols = c("grey", "red")) + 
        ggtitle(paste0("Cluster 15 from all neurons + glia snRNAseq data",
                       "\nin UMAP space (with UL protein)"))+
        labs(color = "Is this cell\nfrom snRNAseq\ncluster 15?"))
dev.off()


# Remove any previous normalization
DefaultAssay(data) <- "RNA"
data_rmUL <- DietSeurat(data,assay = "RNA")

# remove UL proteins from the data
remove_genes <- "ChR2(H134R), mRFP1, UL15, UL18, UL19, UL20, UL21, UL25, YFP"
remove_genes <- unlist(strsplit(remove_genes, split = ","))
remove_genes <- trimws(remove_genes)
data_rmUL <- subset(data_rmUL, 
                       features = setdiff(rownames(data_rmUL), remove_genes))

# SCTransform normalization, PCA and UMAP 
data_rmUL <- SCTransform(data_rmUL, vst.flavor = "v2")
data_rmUL <- RunPCA(data_rmUL,npcs = 20) %>%
  RunHarmony(assay.use="SCT",
             group.by.vars = "library_prep_batch",
             kmeans_init_nstart=20,
             kmeans_init_iter_max=100) %>%
  RunUMAP(reduction = "harmony", dims = 1:20)


pdf("Cluster_15_in_full_data_removeULproteins.pdf",
    width = 12,
    height = 7)
print(DimPlot(data_rmUL, 
              group.by = "clus15", 
              raster = F,
              order = T,
              cols = c("grey", "red")) + 
        ggtitle(paste0("Cluster 15 from all neurons + glia snRNAseq data in",
                       "\nUMAP space (UL proteins were removed from this analysis)"))+
        labs(color = "Is this cell\nfrom snRNAseq\ncluster 15?"))
dev.off()



# END --------------------------------------------------------------------------