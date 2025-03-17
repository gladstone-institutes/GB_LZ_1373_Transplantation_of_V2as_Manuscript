#!/usr/bin/env Rscript
#run locally on Ayushi's laptop

#load required packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)

#set all directory paths
basedir <- paste0("~/Dropbox (Gladstone)/GB-LZ-1373/results/",
                  "09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/",
                  "09_clustering_cell_type_no_CB_no_DF_HumanCells/",
                  "20PC_0.6res_rmBatchEffect_rmCluster16/")
setwd(basedir)

#load the Seurat object
sn_dat_20PC_0.6res <- readRDS("smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_rmBatchEffect_rmCluster16_clustered_and_cell_typed.rds")


# generate heatmap -------------------------------------------------------------
# heatmap for top 10 marker genes for each cluster
de_markers <- read.csv("smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_rmBatchEffect_rmCluster16_FindAllMarkers_marker_genes_per_cluster.csv")
de_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
pdf(file = "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_rmBatchEffect_rmCluster16_heatmap.pdf",
    width = 28,
    height = 18)
print(DoHeatmap(sn_dat_20PC_0.6res, features = top10$gene) + 
        scale_fill_gradientn(colors = c("#33BBEE", "white", "#EE3377")) + 
        ggtitle("Heatmap of marker genes for 20PC_0.6res_rmBatchEffect_rmCluster16 clustering\n(single-nucleus RNA-seq human cells)"))
dev.off()



# dotplots --------------------------------------------------------------------- 
features_df <- read.csv("~/Dropbox (Gladstone)/GB-LZ-1373/assets/Science_Submission_DotPlots_AA_v2.csv")
#The following requested variables were not found: TBFOX3
features_df$gene_symbol <- trimws(features_df$gene_symbol)
features_df <- features_df[features_df$gene_symbol != "TBFOX3",] 
cluster_order <- read.csv("~/Dropbox (Gladstone)/GB-LZ-1373/assets/Science_Submission_DotPlots_AA_cluster_order_v2.csv")
dotplot_id <- unique(features_df$celltype)

# re order clusters as suggested by Lana
sn_dat_20PC_0.6res$seurat_clusters <- factor(sn_dat_20PC_0.6res$seurat_clusters, 
                                             levels = cluster_order$cluster_order)
Idents(sn_dat_20PC_0.6res) <- sn_dat_20PC_0.6res$seurat_clusters

# generate the dot plots
for(i in dotplot_id){
  features_to_plot <- unique(toupper(features_df[features_df$celltype == i, "gene_symbol"]))
  pdf(file = paste0("smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_rmBatchEffect_rmCluster16_",
                    i,"_dotplot.pdf"),
      width = 10)
  print(DotPlot(sn_dat_20PC_0.6res, features = features_to_plot, 
                cols = c("#33BBEE", "#EE3377")) + 
          RotatedAxis() + coord_flip())
  dev.off()
}



# feature plots --------------------------------------------------------------------- 
features_to_plot <- unique(c("SLC17A6","SLC17A8", "GAD1","GAD2","PAX2","SLC6A5", "MAF", "MAFA",
             "PAM", "TAC1", "TAC3", "NMU", "USP9Y", "RBFOX3", "CPAMD8", "SLC6A11", "PAX2", "XRCC5"))
pdf("smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_rmBatchEffect_rmCluster16_featureplots.pdf")
for(i in features_to_plot){
  print(FeaturePlot(sn_dat_20PC_0.6res, 
                    features= i, 
                    raster = FALSE, 
                    order = TRUE, 
                    label = FALSE, 
                    reduction = "umap"))
}
dev.off()



# >0 expression for percentage calculation and violin plot for all marker genes --------
percent_exp <- data.frame(gene_id = character(),
                          percentage_of_cells_with_greater_than_zero_expression = numeric(),
                          stringsAsFactors = F)
for(i in unique(toupper(c(features_df$gene_symbol, features_to_plot, "SLC17A6")))){
  x <- FetchData(sn_dat_20PC_0.6res, vars = i)[,1]
  high_exp_cells_percent <- round(((length(x[x>0])*100)/length(x)),2)
  percent_exp[nrow(percent_exp)+1,] <- c(i,high_exp_cells_percent)
}
write.csv(percent_exp,
          "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_rmBatchEffect_rmCluster16_percentage_expression_of_marker_genes.csv",
          row.names = F)



# cells per sample per cluster -------------------------------------------------
##the number of cells in each mouse in each cluster
all_metadata <- sn_dat_20PC_0.6res[[]]

#make a table of the required data
smp_cluster_counts <- unique(all_metadata %>%
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
colnames(smp_cluster_counts)[1:2] <- c("SampleName","SeuratClusterID")
smp_cluster_counts <- smp_cluster_counts[order(smp_cluster_counts$SampleName),]
write.csv(smp_cluster_counts, 
          file = "counts_per_sample_per_cluster.csv",
          row.names = FALSE)



# cells expressing UL18, UL19 and UL25 for each animal -------------------------
##the number of cells in each mouse in each cluster
all_metadata <- sn_dat_20PC_0.6res[[]]

expr_df <- FetchData(sn_dat_20PC_0.6res, vars = c("UL18", "UL19", "UL25"))

all_metadata <- cbind(all_metadata, expr_df)

# Summarize the number of cells per animalID with UL gene exp > 0
UL_cell_counts_per_sample <- all_metadata %>%
  group_by(sample_name) %>%
  summarize(cells_per_sample = n(),
            UL18_cellcount = sum(UL18 > 0),
            UL19_cellcount = sum(UL19 > 0),
            UL25_cellcount = sum(UL25 > 0),
            any_UL_cellcount = sum(UL25 > 0 | UL18 > 0 | UL19 > 0)) %>%
  as.data.frame()

write.csv(UL_cell_counts_per_sample, 
          file = "cell_count_per_sample_with_UL_gene_exp_greater_than_zero.csv",
          row.names = FALSE)


# END --------------------------------------------------------------------------