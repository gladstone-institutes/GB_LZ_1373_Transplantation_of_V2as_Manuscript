#!/usr/bin/env Rscript
###############################################################################
## Project ID: GB-LZ-1373
## Authors: Ayushi Agrawal
## Script Goal: Merge Seurat objects for sample 8 and single-cell data from 
##          previous experiment; visualize for batch effects.
## Example running this script on wynton:
##    Rscript 08_02_merge_and_visualize_sample_8_old_sample_grch38_no_cellbender_no_doubletfinder.R
###############################################################################

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)

#set all directory paths
basedir <- paste0("/gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-",
                  "scrnaseq-snrnaseq-hg38-may-2023/results/",
                  "08_seurat_analysis_samples_8_and_oldsample_grch38_no_cellbender_no_doubletfinder")
indir <- file.path(basedir,"01_qc")
outdir <- basedir
oldsampledir <- paste0("/gladstone/bioinformatics/projects/lz-1266-lana-zholudeva-deepak-srivastava-scrnaseq",
                       "-hg38-jan-2022/results/05_seurat_analysis_newQCcutoffs_grch38_chr2_yfp_prv/data/01_qc")
setwd(indir)
  
#create the output directory if it doesn't exist
if(!dir.exists(file.path(outdir, "02_merge_and_visualize")))
  dir.create(file.path(outdir, "02_merge_and_visualize"), recursive = T)

#read in all the 26 sample datasets
file_paths <- c(file.path(indir,list.files(pattern = ".rds")),
                file.path(oldsampledir,list.files(path = oldsampledir,pattern = "^LZ-3672_0[1-4]_newQCcutoffs_prv\\.rds$")))
file_names <-  gsub(pattern = "\\.rds$", 
                    replacement = "", 
                    x = basename(file_paths))
data_list <- lapply(file_paths, readRDS)
names(data_list) <- file_names 

#merge all the datasets
data <- merge(data_list[[1]], 
              y=unlist(data_list[2:length(data_list)]), 
              add.cell.ids=file_names,
              project = "LZ-1373-snRNAseq", 
              merge.data=FALSE)
data$library_name <- sub("_[^_]*$", "", colnames(data))
data$data_type <- "single-cell"
data$species <- "human"
data$prv_used <- "no"
data$library_prep_batch <- ifelse(data$sample_name == "LZ4087_08", "new_batch_2", "old_batch")
data$experiment <- ifelse(data$sample_name == "LZ4087_08", "new", "old")
data$library_prep_date <- ifelse(data$sample_name == "LZ4087_08", "new_unknown", "11/03/21")
saveRDS(data, file = paste0(outdir,"/02_merge_and_visualize/merged_data_samples_8_and_oldsample_grch38.rds"))
print("***** Dataset merge completed! *****")


#read in the merged dataset
data <- readRDS(paste0(outdir,"/02_merge_and_visualize/merged_data_samples_8_and_oldsample_grch38.rds"))

#normalization using SCTransform
data_sct <- SCTransform(data, method="glmGamPoi")
saveRDS(data_sct, file = paste0(outdir,"/02_merge_and_visualize/merged_data_samples_8_and_oldsample_grch38_post_sct.rds"))

#Perform and PCA on the SCTransformed data
data_sct_processed <- RunPCA(object = data_sct, assay = "SCT")
pdf(file = file.path(outdir, "/02_merge_and_visualize/merged_data_samples_8_and_oldsample_grch38_PCA_plots.pdf"),
    width = 12,
    height = 7)
print(ElbowPlot(data_sct_processed, ndims = 50))
print(ElbowPlot(data_sct_processed, ndims = 20))
print(VizDimLoadings(data_sct_processed, dims = 1:2, reduction = "pca"))
print(DimPlot(data_sct_processed, reduction = "pca", raster = FALSE))
print(DimHeatmap(data_sct_processed, dims = 1:2, cells = 500, balanced = TRUE))
print(DimHeatmap(data_sct_processed, dims = 3:4, cells = 500, balanced = TRUE))
dev.off()

#Run UMAP on the data
data_sct_processed <- data_sct_processed %>% RunUMAP(., dims = 1:10, verbose = TRUE)
saveRDS(data_sct_processed, file = paste0(outdir,"/02_merge_and_visualize/",
                                          "merged_data_samples_8_and_oldsample_grch38_post_sct_processed.rds"))

#Generate visualizations using the various metadata
plot_metadata <- function(sc_dat, metadata_variable, pdf_height){
  DimPlot(sc_dat, 
          raster = FALSE, 
          order = TRUE, 
          label = FALSE, 
          group.by=metadata_variable, 
          reduction = "umap") %>%
    ggsave(file = file.path(outdir, paste0("02_merge_and_visualize/merged_data_samples_8_and_oldsample_grch38_", 
                                           metadata_variable, "_umap.pdf")),
           plot = .,
           width = 12,
           height = 7)
  DimPlot(sc_dat, 
          raster = FALSE, 
          order = TRUE, 
          label = FALSE, 
          group.by=metadata_variable,
          split.by=metadata_variable, 
          reduction = "umap",
          ncol = 2) %>%
    ggsave(file = file.path(outdir, paste0("02_merge_and_visualize/merged_data_samples_8_and_oldsample_grch38_", 
                                           metadata_variable,"_split_umap.pdf")),
           plot = .,
           width = 12,
           height = pdf_height)
}

plot_metadata(data_sct_processed, "sample_name", 12)
plot_metadata(data_sct_processed, "library_name", 12)
plot_metadata(data_sct_processed, "library_prep_batch", 7)
plot_metadata(data_sct_processed, "experiment", 7)
plot_metadata(data_sct_processed, "library_prep_date", 7)


FeaturePlot(data_sct_processed, 
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            features="nCount_RNA", 
            reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("02_merge_and_visualize/merged_data_samples_8_and_oldsample_grch38", 
                                 "_nCount_RNA_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
FeaturePlot(data_sct_processed, 
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            features="nFeature_RNA", 
            reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("02_merge_and_visualize/merged_data_samples_8_and_oldsample_grch38", 
                                 "_nFeature_RNA_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
FeaturePlot(data_sct_processed, 
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            features="percent.mt", 
            reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("02_merge_and_visualize/merged_data_samples_8_and_oldsample_grch38", 
                                 "_percent.mt_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
FeaturePlot(data_sct_processed, 
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            features="nCount_SCT", 
            reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("02_merge_and_visualize/merged_data_samples_8_and_oldsample_grch38", 
                                 "_nCount_SCT_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
FeaturePlot(data_sct_processed, 
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            features="nFeature_SCT", 
            reduction = "umap") %>%
  ggsave(file = file.path(outdir, paste0("02_merge_and_visualize/merged_data_samples_8_and_oldsample_grch38", 
                                         "_nFeature_SCT_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)
FeaturePlot(data_sct_processed, 
            raster = FALSE, 
            order = TRUE, 
            label = FALSE, 
            features="MBP", 
            reduction = "umap") %>%
  ggsave(file = file.path(outdir, paste0("02_merge_and_visualize/merged_data_samples_8_and_oldsample_grch38", 
                                         "_MBPgene_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)

#save the session info
writeLines(capture.output(sessionInfo()), 
           file.path(outdir, "sessionInfo_step02_merge_and_visualize_sample_8_and_oldsample_grch38.txt"))


print("***** Merge, normalization and visualization completed! *****")


########################## END ##########################  

