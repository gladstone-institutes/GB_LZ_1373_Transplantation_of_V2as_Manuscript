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
                  "20PC_0.6res_rmBatchEffect/")
setwd(basedir)

#load the Seurat object
sn_dat_20PC_0.6res <- readRDS("smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_rmBatchEffect_clustered_and_cell_typed.rds")



# generate heatmap -------------------------------------------------------------
# heatmap for top 10 marker genes for each cluster
de_markers <- read.csv("smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_rmBatchEffect_FindAllMarkers_marker_genes_per_cluster.csv")
de_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
pdf(file = "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_rmBatchEffect_heatmap.pdf",
    width = 28,
    height = 18)
print(DoHeatmap(sn_dat_20PC_0.6res, features = top10$gene) + 
        scale_fill_gradientn(colors = c("#33BBEE", "white", "#EE3377")) + 
        ggtitle("Heatmap of marker genes for 20PC_0.6res_rmBatchEffect clustering\n(single-nucleus RNA-seq human cells)"))
dev.off()



# dotplots --------------------------------------------------------------------- 
# general population features
general_features <- c("tubb3",	"rbfox3",	"map2", "sox2",	"pax6",	"MKI67", "pax3",	"pax7", "Nkx6-1",	"Prdm8",
                      "slc1a3",	"SPARC", "ATP1B2", "GFAP",	"S100B",	"AQP4","CPAMD8",	"CCDC85A",	
                      "HAP1",	"GABBR2",	"SLC6A11", "CSPG4",	"PLP1",	"PDGFRA", "MOG","MAG", "mpz","pmp22",
                      "DLX6", "OLIG3", "XRCC5")
pdf(file = "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_rmBatchEffect_general_population_dotplot.pdf",
    width = 10)
print(DotPlot(sn_dat_20PC_0.6res, features = toupper(general_features), 
              cols = c("#33BBEE", "#EE3377")) + 
        RotatedAxis() + coord_flip())
dev.off()

# dorsal vs ventral features
dorsa_ventral_features <- c("Lbx1","Tlx3","foxd3","Satb1",	"Tlx3","Lbx1","Tlx3","Lbx1","EVX1",
                            "Lbx1","Satb1","PAX2", "Lhx5","En1","Foxp2", "foxd3","VSX1", "DLL3", "DLL4",
                            "Vsx2","Lhx3","sox14","sox21", "ONECUT2","Gata2","ISL1", "MNX1", "PCP4", "MNX1")
#could not find genes: "PHOX2B", "CALRETININ","CALBINDIN", "BHLHB5", "PHOX2A", , "DBX1", "ISL2"
pdf(file = "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_rmBatchEffect_dorsal_ventral_population_dotplot.pdf",
    width = 10)
print(DotPlot(sn_dat_20PC_0.6res, features = unique(toupper(dorsa_ventral_features)), 
              cols = c("#33BBEE", "#EE3377")) +  
        RotatedAxis() + coord_flip())
dev.off()

# Excitatory vs Inhibitory vs Peptidergic
ex_in_pep_features <- c("SLC17A7",	"SLC17A6",	"SLC17A8",	"GAD1",	"GAD2",	"PAX2",	"SLC6A5", "ADARB2",	"PDYN",	"RORB",	
                        "MAFA",	"TAC1",	"SOX5",	"PDE11A", "ZFHX3",	"NFIX", "MAF",	"MAFA", "PAM",	"TAC1",	"TAC3",	"NMU")
pdf(file = "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_rmBatchEffect_excitatory_inhibitory_peptidergic_dotplot.pdf",
    width = 10)
print(DotPlot(sn_dat_20PC_0.6res, features = unique(toupper(ex_in_pep_features)), 
              cols = c("#33BBEE", "#EE3377")) +  
        RotatedAxis() + coord_flip())
dev.off()

# HOX genes
hox_genes_suppl_fig <- c("LMX1A","FOXA2","SIM1","PAX5","GBX2","HOXA1","HOXA2","HOXA3","HOXB2","HOXB3",
                         "HOXD1","HOXD3","HOXA4","HOXA5",'HOXB4',"HOXB5","HOXC4","HOXC5","HOXD4","HOXA6","HOXA7","HOXB6",
                         "HOXB7","HOXC6","HOXB8","HOXC8","HOXD8","HOXA9","HOXB9","HOXC9","HOXA10",
                         "HOXB13")
#could not find genes: "FOXG1", "HOXD9","HOXC10","HOXB1", "HOXA8", "HOXA11", "HOXC11", 
#"HOXD13", "HOXA12", "HOXC12", "HOXD12", "HOXD10", "HOXA13", "HOXC13"
pdf(file = "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_rmBatchEffect_hox_genes_dotplot.pdf",
    width = 10)
print(DotPlot(sn_dat_20PC_0.6res, features = unique(toupper(hox_genes_suppl_fig)), 
              cols = c("#33BBEE", "#EE3377")) +  
        RotatedAxis() + coord_flip())
dev.off()



# feature plots --------------------------------------------------------------------- 
features_to_plot <- unique(c("SLC17A7","SLC17A6","SLC17A8", "GAD1","GAD2","PAX2","SLC6A5", "MAF", "MAFA",
             "PAM", "TAC1", "TAC3", "NMU", "USP9Y", "RBFOX3", "CPAMD8", "SLC6A11", "PAX2", "XRCC5"))
pdf("smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_rmBatchEffect_featureplots.pdf")
for(i in features_to_plot){
  print(FeaturePlot(sn_dat_20PC_0.6res, 
                    features= i, 
                    raster = FALSE, 
                    order = TRUE, 
                    label = FALSE, 
                    reduction = "umap"))
}
dev.off()


##################
# >0 expression for percentage calculation and violin plot for all marker genes
##################
percent_exp <- data.frame(gene_id = character(),
                          percentage_of_cells_with_greater_than_zero_expression = numeric(),
                          stringsAsFactors = F)
for(i in unique(toupper(c(general_features, dorsa_ventral_features, ex_in_pep_features, 
                          features_to_plot, hox_genes_suppl_fig, "SLC17A6")))){
  x <- FetchData(sn_dat_20PC_0.6res, vars = i)[,1]
  high_exp_cells_percent <- round(((length(x[x>0])*100)/length(x)),2)
  percent_exp[nrow(percent_exp)+1,] <- c(i,high_exp_cells_percent)
}
write.csv(percent_exp,
          "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_rmBatchEffect_percentage_expression_of_marker_genes.csv",
          row.names = F)



#################### END ####################
