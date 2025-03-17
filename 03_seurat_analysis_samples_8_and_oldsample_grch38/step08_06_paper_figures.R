#!/usr/bin/env Rscript
#run locally on Ayushi's laptop

#load required packages --------------------------------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)


#set all directory paths
basedir <- paste0("~/Dropbox (Gladstone)/GB-LZ-1373/results/",
                  "08_seurat_analysis_samples_8_and_oldsample_grch38_no_cellbender_no_doubletfinder")



# define the function ----------------------------------------------------------
generate_paper_figure2 <- function(dat, res){
  
  #Fig 2B: featureplot ---------------------------------------------------------
  features_to_plot <- c("VSX1","VSX2")
  pdf(paste0("samples_8_and_oldsample_grch38_singlecell_15PC_",res,
             "res_featureplot_fig2b.pdf"),
      width = 15)
  print(FeaturePlot(dat,
                    features = features_to_plot,
                    raster = FALSE,
                    order = TRUE,
                    label = FALSE,
                    reduction = "umap",
                    ncol = 2))
  dev.off()
  
  
  #Fig 2C: heatmap -------------------------------------------------------------
  #heatmap for top 10 marker genes for each cluster
  de_markers <- read.csv("samples_8_and_oldsample_grch38_FindAllMarkers_marker_genes_per_cluster.csv")
  de_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
  pdf(file = paste0("samples_8_and_oldsample_grch38_singlecell_15PC_",
                    res,"res_heatmap_fig2c.pdf"),
      width = 12,
      height = 15)
  print(DoHeatmap(dat, features = top10$gene) +
          scale_fill_gradientn(colors = c("#33BBEE", "white", "#EE3377")))
  dev.off()
  
  
  #Fig 2D-F: dot plots -----------------------------------------------------------
  #HOX genes rostro-caudal phenotype
  hox_genes_rostro_caudal <- c('FOXG1', 'LMX1A', 'FOXA2', 'SIM1', 'HOXA1', 'HOXA2', 'HOXA3', 'HOXA4',
                               'HOXA5', 'HOXB4', 'HOXB6', 'HOXB7', 'HOXC6', 'HOXC8', 'HOXD8', 'HOXC9',
                               'HOXA9', 'HOXD9', 'HOXC10', 'HOXD10', 'HOXA13', 'HOXB13', 'HOXC13')
  couldnot_find_genes <- c("HOXA13", "HOXB13", "HOXC10", "HOXD10", "HOXC13")
  hox_genes_rostro_caudal <- hox_genes_rostro_caudal[!(hox_genes_rostro_caudal %in% couldnot_find_genes)]
  pdf(file = paste0("samples_8_and_oldsample_grch38_singlecell_15PC_", res,
                    "res_hox_genes_rostro-caudal_phenotype_dotplot_scaled_fig2d.pdf"))
  print(DotPlot(dat, 
                features = unique(toupper(hox_genes_rostro_caudal)),
                cols = c("#33BBEE", "#EE3377")) +
          RotatedAxis() + coord_flip())
  dev.off()
  
  #general population features
  general_features <- c("TUBB3","rbfox3",	"map2", "sox2",	"pax6",	"MKI67", "pax3", "pax7",
                        "Nkx6-1",	"Prdm8", "FGFR3",	"slc1a3", "GFAP",	"S100B",	"AQP4",
                        "CSPG4",	"PLP1",	"PDGFRA","Olig3",	"MOG",	"MAG", "mpz",	"pmp22"	)
  pdf(file = paste0("samples_8_and_oldsample_grch38_singlecell_15PC_", res,
                    "res_general_population_dotplot_scaled_fig2e.pdf"))
  print(DotPlot(dat, 
                features = toupper(general_features),
                cols = c("#33BBEE", "#EE3377")) +
          RotatedAxis() + coord_flip())
  dev.off()
  
  #dorsal vs ventral features
  dorsa_ventral_features <- c("Lbx1",	"Tlx3",	"Satb1", "PHOX2A",	"PAX2",	"PHOX2B",
                              "Evx1",	"Lhx5",	"En1", "Foxp2",	"Foxp1",
                              "VSX1",	"DLL3",	"DLL4", "Vsx2",	"Lhx3",	"sox14",
                              "sox21", "ONECUT2", "Gata2","ISL1", "MNX1",	"ISL2",
                              "PCP4", "MNX1")
  #could not find genes: "Calretinin","calbindin", "Bhlhb5", "Ptx2", "HB9"
  pdf(file = paste0("samples_8_and_oldsample_grch38_singlecell_15PC_", res,
                    "res_dorsal_ventral_population_dotplot_scaled_fig2f.pdf"))
  print(DotPlot(dat, 
                features = unique(toupper(dorsa_ventral_features)),
                cols = c("#33BBEE", "#EE3377")) +
          RotatedAxis() + coord_flip())
  dev.off()
  
  
  # session information
  writeLines(
    capture.output(sessionInfo()),
    file.path(outdir, "sessionInfo_AA_local.txt")
  )
  
}




# run for 15 PC and 0.08 res ----------------------------------------------------
outdir <- file.path(basedir, "04_clustering_cell_type/15PC_0.08res")

setwd(outdir)

# load the Seurat object
sc_dat_15PC_0.08res <- readRDS("samples_8_and_oldsample_grch38_clustered_and_cell_typed.rds")

# generate all figures
generate_paper_figure2(sc_dat_15PC_0.08res, 0.08)




# run for 15 PC and 0.1 res ----------------------------------------------------
outdir <- file.path(basedir, "04_clustering_cell_type/15PC_0.1res")

setwd(outdir)

# load the Seurat object
sc_dat_15PC_0.1res <- readRDS("samples_8_and_oldsample_grch38_clustered_and_cell_typed.rds")

# generate all figures
generate_paper_figure2(sc_dat_15PC_0.1res, 0.1)



# END -------------------------------------------------------------------------