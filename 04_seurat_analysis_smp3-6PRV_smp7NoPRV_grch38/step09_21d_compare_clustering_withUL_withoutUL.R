# Setup ------------------------------------------------------------------------
# load libraries
library(Seurat)
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)


# set working directory
setwd("~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF")


# Inputs -----------------------------------------------------------------------
dat_UL <- readRDS("09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.6res_rmBatchEffect/smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_rmBatchEffect_clustered_and_cell_typed.rds")
dat_UL_rmclus16 <- readRDS("09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.6res_rmBatchEffect_rmCluster16/smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_rmBatchEffect_rmCluster16_clustered_and_cell_typed.rds")

dat_noUL <- readRDS("15_clustering_cell_type_no_CB_no_DF_HumanCells_NoPRV/20PC_0.6res_rmBatchEffect/smp3-7_grch38_noCB_noDF_HumanCells_rmBatchEffect_NoPRV_20PC_0.6res_clustered_and_cell_typed.rds")
dat_noUL_rmclus15 <- readRDS("15_clustering_cell_type_no_CB_no_DF_HumanCells_NoPRV/20PC_0.6res_rmBatchEffect_rmCluster15/smp3-7_grch38_noCB_noDF_HumanCells_rmBatchEffect_NoPRV_20PC_0.6res_rmCluster15_clustered_and_cell_typed.rds")


# Main -------------------------------------------------------------------------
# before removing artifact cluster
# extract clusters
cluster_UL <- dat_UL@meta.data$seurat_clusters
names(cluster_UL) <- rownames(dat_UL@meta.data)
cluster_noUL <- dat_noUL@meta.data$seurat_clusters
names(cluster_noUL) <- rownames(dat_noUL@meta.data)

# create a data frame with both clustering results
cluster_noUL <- cluster_noUL[names(cluster_UL)]

# create a data frame with both clustering results
cluster_data <- data.frame(cluster_UL = cluster_UL, 
                           cluster_noUL = cluster_noUL)

# Transform the matrix in long format
data = cluster_data %>% group_by(cluster_UL, cluster_noUL) %>% 
  summarise(n=log10(n())) %>% as.data.frame()

# plot the data
pdf("15_clustering_cell_type_no_CB_no_DF_HumanCells_NoPRV/cluster_membership_20PC_0.6res_rmBatchEffect_withUL_vs_withoutUL.pdf",
    width = 22,
    height = 14)
print(ggplot(data, aes(x = cluster_UL, y = cluster_noUL, fill = n)) +
        geom_tile(color = "black") +
        geom_text(aes(label = round(10^n,2)), color = "black", size = 4) +
        scale_fill_gradient(low = "#e0f3f8", high = "red") +
        labs(x = "Clusters at 0.6 resolution (with UL transcripts)", 
             y = "Clusters at 0.6 resolution (without UL transcripts)", 
             fill = "log10(number of cells)", 
             title = "Number of cells in the same cluster at 0.6 resolution (with/without UL transcripts)") +
        theme_bw())
dev.off()




# after removing artifact cluster
# extract clusters
cluster_UL_rmclus <- dat_UL_rmclus16@meta.data$seurat_clusters
names(cluster_UL_rmclus) <- rownames(dat_UL_rmclus16@meta.data)
cluster_noUL_rmclus <- dat_noUL_rmclus15@meta.data$seurat_clusters
names(cluster_noUL_rmclus) <- rownames(dat_noUL_rmclus15@meta.data)

# create a data frame with both clustering results
cluster_UL_rmclus <- cluster_UL_rmclus[names(cluster_noUL_rmclus)]

# create a data frame with both clustering results
cluster_data_rmclus <- data.frame(cluster_UL_rmclus = cluster_UL_rmclus, 
                                  cluster_noUL_rmclus = cluster_noUL_rmclus)

# Transform the matrix in long format
data_rmclus = cluster_data_rmclus %>% 
  group_by(cluster_UL_rmclus, cluster_noUL_rmclus) %>% 
  summarise(n=log10(n())) %>% as.data.frame()

# plot the data
pdf("15_clustering_cell_type_no_CB_no_DF_HumanCells_NoPRV/cluster_membership_20PC_0.6res_rmBatchEffect_rmclus_withUL_vs_withoutUL.pdf",
    width = 22,
    height = 14)
print(ggplot(data_rmclus, aes(x = cluster_UL_rmclus, y = cluster_noUL_rmclus, fill = n)) +
        geom_tile(color = "black") +
        geom_text(aes(label = round(10^n,2)), color = "black", size = 4) +
        scale_fill_gradient(low = "#e0f3f8", high = "red") +
        labs(x = "Clusters at 0.6 resolution (with UL transcripts; after removing cluster 16)", 
             y = "Clusters at 0.6 resolution (without UL transcripts; after removing cluster 15)", 
             fill = "log10(number of cells)", 
             title = "Number of cells in the same cluster at 0.6 resolution after removing low-quality cluster (with/without UL transcripts)") +
        theme_bw())
dev.off()



# Check clusters for sample-driven clusters ---------------------
smp_cluster_counts <- unique(dat_noUL_rmclus15@meta.data %>%
                                      group_by(sample_name) %>%
                                      mutate(total_numbers_of_cells_per_sample = 
                                               n()) %>%
                                      group_by(seurat_clusters, .add=TRUE) %>%
                                      mutate(number_of_cells_per_sample_in_cluster = 
                                               n()) %>%
                                      select(sample_name,
                                             seurat_clusters,
                                             total_numbers_of_cells_per_sample,
                                             number_of_cells_per_sample_in_cluster))
colnames(smp_cluster_counts)[1:2] <- c("sample_id","cluster_id")
smp_cluster_counts <- smp_cluster_counts[order(smp_cluster_counts$cluster_id),]

smp_cluster_counts <- smp_cluster_counts %>% 
  group_by(cluster_id) %>% 
  mutate(total_numbers_of_cells_per_cluster=sum(number_of_cells_per_sample_in_cluster)) %>% 
  group_by(sample_id, .add = TRUE) %>% 
  mutate(percent_per_sample_in_cluster = round((number_of_cells_per_sample_in_cluster/total_numbers_of_cells_per_cluster), 2), 
         percent_per_sample_in_cluster_abs = (number_of_cells_per_sample_in_cluster/total_numbers_of_cells_per_cluster))

(ggplot(smp_cluster_counts, aes(x = factor(cluster_id), y = percent_per_sample_in_cluster_abs, 
                                        fill = factor(sample_id))) + 
    geom_bar(stat = "identity", position = "stack", color = "lightgrey", size = 0.08) +
    labs(x = "Cluster", y = "Percent", fill = "Sample", 
         title = "Stacked Bar Plot - percentage of cells in a cluster from a sample (0.6 resolution; rm clus 15; no UL)") +
    theme_bw() +
    theme(text = element_text(size = 20))) %>% 
  ggsave("15_clustering_cell_type_no_CB_no_DF_HumanCells_NoPRV/20PC_0.6res_rmBatchEffect_rmCluster15/percent_cells_in_cluster_from_sample.pdf",.,
         width = 16)

(ggplot(smp_cluster_counts, aes(x = factor(cluster_id), y = number_of_cells_per_sample_in_cluster, 
                                       fill = factor(sample_id))) + 
    geom_bar(stat = "identity", position = "stack", color = "lightgrey", size = 0.08) +
    labs(x = "Cluster", y = "Count", fill = "Sample", 
         title = "Stacked Bar Plot - number of cells per sample in a cluster (0.6 resolution; rm clus 15; no UL)") +
    theme_bw() +
    theme(text = element_text(size = 20))) %>%
  ggsave("15_clustering_cell_type_no_CB_no_DF_HumanCells_NoPRV/20PC_0.6res_rmBatchEffect_rmCluster15/cells_per_sample_barplot.pdf",.,
         width = 12)


# END --------------------------------------------------------------------------