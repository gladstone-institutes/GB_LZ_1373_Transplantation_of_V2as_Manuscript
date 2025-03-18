library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)


basedir <- "~/Dropbox (Gladstone)/GB-LZ-1373/results/"
setwd(paste0(basedir, 
             "10_subcluster_neurons_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/",
             "02_neuron_subclustering/30PC_0.10res"))

neurons <- readRDS("neurons_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_30PC_0.10res_clustered_and_cell_typed.rds")
all_data <- readRDS(paste0(basedir, 
                           "09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/",
                           "09_clustering_cell_type_no_CB_no_DF_HumanCells/",
                           "20PC_0.6res_rmBatchEffect_rmCluster16/",
                           "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_",
                           "rmBatchEffect_rmCluster16_clustered_and_cell_typed.rds"))



# 1. Where is nueuron subcluster 5 in all data? ---------------------------------------------------------
cell_ids_subclus5 <- WhichCells(neurons, idents = 5)
all_data$neuron_0.1_subclus5 <- ifelse(rownames(all_data@meta.data) %in% cell_ids_subclus5, "Yes", "No")
table(all_data$neuron_0.1_subclus5)
pdf("neuron_subcluster_5_30PC_0.1res_in_all_snRNAseqData.pdf",
    width = 12,
    height = 7)
print(DimPlot(all_data, 
              group.by = "neuron_0.1_subclus5", 
              raster = F,
              order = T,
              cols = c("grey", "red")) + 
        ggtitle("Neuron subcluster 5 (30PC_0.1res) cells \nin all snRNAseqData UMAP space")+
        labs(color = "Is this cell\nin neuron\nsubcluster 5?"))
dev.off()




# 2. Where is nueuron subcluster 5 in all data? ---------------------------------------------------------
cell_ids_subclus2 <- WhichCells(neurons, idents = 2)
all_data$neuron_0.1_subclus2 <- ifelse(rownames(all_data@meta.data) %in% cell_ids_subclus2, "Yes", "No")
table(all_data$neuron_0.1_subclus2)
pdf("neuron_subcluster_2_30PC_0.1res_in_all_snRNAseqData.pdf",
    width = 12,
    height = 7)
print(DimPlot(all_data, 
              group.by = "neuron_0.1_subclus2", 
              raster = F,
              order = T,
              cols = c("grey", "red")) + 
        ggtitle("Neuron subcluster 2 (30PC_0.1res) cells \nin all snRNAseqData UMAP space")+
        labs(color = "Is this cell\nin neuron\nsubcluster 2?"))
dev.off()




# 3. Take out the UL (UL18, 19, 25 and mRFP) transcripts from the sub clustered 
#    data set and see if this “Cluster 15” remains a unique cluster, or if the 
#    cells are now dispersed throughout the dataset -----------------------------------------------------
# visualize cluster 15
neurons$clus15 <- ifelse(neurons$SCT_snn_res.0.6 == 15, "Yes", "No")
pdf("Cluster_15_in_neuron_subcluster_withULproteins.pdf",
    width = 12,
    height = 7)
print(DimPlot(neurons, 
              group.by = "clus15", 
              raster = F,
              order = T,
              cols = c("grey", "red")) + 
        ggtitle(paste0("Cluster 15 from all snRNAseq data in neuron subcluster UMAP space",
                       "\n(with UL protein)"))+
        labs(color = "Is this cell\nfrom snRNAseq\ncluster 15?"))
dev.off()


# Remove any previous normalization
DefaultAssay(neurons) <- "RNA"
neurons_rmUL <- DietSeurat(neurons,assay = "RNA")

remove_genes <- "ChR2(H134R), mRFP1, UL15, UL18, UL19, UL20, UL21, UL25, YFP"
remove_genes <- unlist(strsplit(remove_genes, split = ","))
remove_genes <- trimws(remove_genes)

neurons_rmUL <- subset(neurons_rmUL, 
                       features = setdiff(rownames(neurons_rmUL), remove_genes))

# SCTransform normalization, PCA and UMAP 
neurons_rmUL <- SCTransform(neurons_rmUL, vst.flavor = "v2")
neurons_rmUL <- RunPCA(neurons_rmUL,npcs = 30) %>%
  RunHarmony(assay.use="SCT",
             group.by.vars = "library_prep_batch",
             kmeans_init_nstart=20,
             kmeans_init_iter_max=100) %>%
  RunUMAP(reduction = "harmony", dims = 1:30)


pdf("Cluster_15_in_neuron_subcluster_removeULproteins.pdf",
    width = 12,
    height = 7)
print(DimPlot(neurons_rmUL, 
              group.by = "clus15", 
              raster = F,
              order = T,
              cols = c("grey", "red")) + 
        ggtitle(paste0("Cluster 15 from all snRNAseq data in neuron subcluster UMAP space",
                       "\n(UL proteins were removed from this analysis)"))+
        labs(color = "Is this cell\nfrom snRNAseq\ncluster 15?"))
dev.off()



# END --------------------------------------------------------------------------