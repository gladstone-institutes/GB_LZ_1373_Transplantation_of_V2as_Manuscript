#!/usr/bin/env Rscript
#run locally on Ayushi's laptop

#load required packages
library(Seurat)
library(ggplot2)
library(dplyr)


# full dataset (neurons + glia) snRNA plots ------------------------------------
#set all directory paths
basedir <- file.path("~/Dropbox (Gladstone)/GB-LZ-1373/results",
                     "09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF",
                     "15_clustering_cell_type_no_CB_no_DF_HumanCells_NoPRV",
                     "20PC_0.6res_rmBatchEffect_rmCluster15/")
setwd(basedir)

#load the Seurat object
sn_dat_20PC_0.6res <- readRDS("smp3-7_grch38_noCB_noDF_HumanCells_rmBatchEffect_NoPRV_20PC_0.6res_rmCluster15_clustered_and_cell_typed.rds")

# featureplot
features_to_plot <- toupper(unique(c("rbfox3", 
                                     "GFAP", 
                                     "slc17a6", 
                                     "GAD1", 
                                     "GAD2", 
                                     "usp9y")))
pdf("smp3-7_grch38_noCB_noDF_HumanCells_rmBatchEffect_NoPRV_20PC_0.6res_rmCluster15_featureplots_paper_fig.pdf")
for(i in features_to_plot){
  print(FeaturePlot(sn_dat_20PC_0.6res, 
                    features= i, 
                    raster = FALSE, 
                    order = TRUE, 
                    label = FALSE, 
                    reduction = "umap"))
}
dev.off()

# cells per sample per cluster
smp_cluster_counts <- unique(sn_dat_20PC_0.6res@meta.data %>%
                               group_by(sample_name) %>%
                               mutate(total_numbers_of_cells_per_sample = 
                                        n()) %>%
                               group_by(seurat_clusters, .add=TRUE) %>%
                               mutate(number_of_cells_per_sample_in_cluster = 
                                        n())%>%
                               select(sample_name,
                                      seurat_clusters,
                                      total_numbers_of_cells_per_sample,
                                      number_of_cells_per_sample_in_cluster))
colnames(smp_cluster_counts)[1:2] <- c("SampleName","ClusterID")
smp_cluster_counts <- smp_cluster_counts[order(smp_cluster_counts$ClusterID),]
write.csv(smp_cluster_counts, 
          file = "smp3-7_grch38_noCB_noDF_HumanCells_rmBatchEffect_NoPRV_20PC_0.6res_rmCluster15_counts_per_sample.csv",
          row.names = FALSE)



# neurons only snRNA plots ------------------------------------
#set all directory paths
basedir <- file.path("~/Dropbox (Gladstone)/GB-LZ-1373/results",
                  "10_subcluster_neurons_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF",
                  "05_neuron_subclustering_NoPRV",
                  "30PC_0.30res_rmBatchEffect")
setwd(basedir)

#load the Seurat object
neurons <- readRDS("neurons_smp3-7_grch38_noCB_noDF_NoPRV_30PC_0.30res_rmBatchEffect_clustered_and_cell_typed.rds")

# featureplot
features_to_plot <- toupper(unique(c("rbfox3", 
                                     "slc17a6",
                                     "GAD1",
                                     "GAD2",
                                     "usp9y",
                                     "tph2",
                                     "chat")))
pdf("neurons_smp3-7_grch38_noCB_noDF_NoPRV_30PC_0.30res_rmBatchEffect_featureplots_paper_fig.pdf")
for(i in features_to_plot){
  print(FeaturePlot(neurons, 
                    features= i, 
                    raster = FALSE, 
                    order = TRUE, 
                    label = FALSE, 
                    reduction = "umap"))
}
dev.off()

# cells per sample per cluster
smp_cluster_counts <- unique(neurons@meta.data %>%
                               group_by(sample_name) %>%
                               mutate(total_numbers_of_cells_per_sample = 
                                        n()) %>%
                               group_by(seurat_clusters, .add=TRUE) %>%
                               mutate(number_of_cells_per_sample_in_cluster = 
                                        n())%>%
                               select(sample_name,
                                      seurat_clusters,
                                      total_numbers_of_cells_per_sample,
                                      number_of_cells_per_sample_in_cluster))
colnames(smp_cluster_counts)[1:2] <- c("SampleName","ClusterID")
smp_cluster_counts <- smp_cluster_counts[order(smp_cluster_counts$ClusterID),]
write.csv(smp_cluster_counts, 
          file = "neurons_smp3-7_grch38_noCB_noDF_NoPRV_30PC_0.30res_rmBatchEffect_counts_per_sample.csv",
          row.names = FALSE)


# END --------------------------------------------------------------------------