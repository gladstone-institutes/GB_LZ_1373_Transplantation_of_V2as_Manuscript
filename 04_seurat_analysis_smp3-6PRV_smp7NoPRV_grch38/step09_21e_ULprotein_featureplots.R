# load required libraries
library(Seurat)

# load the dataset
dat <- readRDS(paste0("~/Dropbox (Gladstone)/GB-LZ-1373/results/",
                      "09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/",
                      "15_clustering_cell_type_no_CB_no_DF_HumanCells_NoPRV/",
                      "20PC_0.6res_rmBatchEffect_rmCluster15/",
                      "smp3-7_grch38_noCB_noDF_HumanCells_rmBatchEffect_NoPRV_",
                      "20PC_0.6res_rmCluster15_clustered_and_cell_typed.rds"))


# > colnames(dat@meta.data)
# [1] "orig.ident"                     "nCount_RNA"                     "nFeature_RNA"                   "percent.mt"                    
# [5] "percent.ribo"                   "sample_name"                    "library_name"                   "data_type"                     
# [9] "species"                        "prv_used"                       "library_prep_batch"             "experiment"                    
# [13] "nCount_sctwithultranscripts"    "nFeature_sctwithultranscripts"  "SCT_norm_expr_ChR2.H134R."      "SCT_norm_expr_mRFP1"           
# [17] "SCT_norm_expr_UL25"             "SCT_norm_expr_UL21"             "SCT_norm_expr_UL19"             "SCT_norm_expr_UL18"            
# [21] "SCT_norm_expr_UL15"             "SCT_norm_expr_YFP"              "SCT_norm_expr_UL20"             "nCount_SCT"                    
# [25] "nFeature_SCT"                   "SCT_snn_res.0.6"                "seurat_clusters"                "cell_type_sctype"              
# [29] "seurat_cluster_celltype_sctype"


# visualize data
pdf(paste0("~/Dropbox (Gladstone)/GB-LZ-1373/results/",
           "09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/",
           "15_clustering_cell_type_no_CB_no_DF_HumanCells_NoPRV/",
           "20PC_0.6res_rmBatchEffect_rmCluster15/",
           "smp3-7_grch38_noCB_noDF_HumanCells_rmBatchEffect_NoPRV_",
           "20PC_0.6res_rmCluster15_UL_proteins_featureplot.pdf"))
for (meta in c("SCT_norm_expr_ChR2.H134R.",
               "SCT_norm_expr_mRFP1",
               "SCT_norm_expr_YFP",
               "SCT_norm_expr_UL15",
               "SCT_norm_expr_UL18",
               "SCT_norm_expr_UL19",
               "SCT_norm_expr_UL20",
               "SCT_norm_expr_UL21",
               "SCT_norm_expr_UL25")) {
  print(FeaturePlot(dat, 
                    features= meta, 
                    raster = FALSE, 
                    order = TRUE, 
                    label = TRUE, 
                    reduction = "umap"))
}
dev.off()


## END 