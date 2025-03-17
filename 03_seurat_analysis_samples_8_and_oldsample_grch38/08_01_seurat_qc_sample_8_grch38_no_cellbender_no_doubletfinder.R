#!/usr/bin/env Rscript
###############################################################################
## Project ID: GB-LZ-1373
## Authors: Ayushi Agrawal
## Script Goal: Create Seurat object from cell ranger count filtered matrices 
##              samples 3,4,5,6, and 7; perform QC on both samples
## Example running this script on wynton:
##    Rscript 08_01_seurat_qc_sample_8_grch38_no_cellbender_no_doubletfinder.R
###############################################################################

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)

#set all directory paths
basedir <- paste0("/gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-",
                  "deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/results")
indir <- paste0(basedir,"/01_cellranger_count/")
outdir <- paste0(basedir,"/08_seurat_analysis_samples_8_and_oldsample_grch38_no_cellbender_no_doubletfinder")
path_postfix_filtered_matrix <- "/outs/filtered_feature_bc_matrix"

#create the output directory
if(!dir.exists(file.path(outdir, "01_qc")))
  dir.create(file.path(outdir, "01_qc"), recursive = TRUE)
setwd(outdir)


#get the list of all rat/human spinal cord samples
samples_list <- c("LZ4087_08A_human","LZ4087_08B_human")

#create a data frame to record cell count before and after filtering per sample
cells_per_sample <- data.frame(sample=character(), 
                               filter_stage=character(), 
                               number_of_cells=numeric(),
                               nFeature_RNA_min_cutoff=numeric(),
                               nFeature_RNA_max_cutoff=numeric(),
                               percent.mt_cutoff=numeric(), 
                               stringsAsFactors = FALSE)

#read in data for each sample and perform QC
for(i in 1:length(samples_list)){
  print(paste0("Analyzing sample: ", samples_list[i], "----------------------"))
  obj.name <- samples_list[i]
  datadir <- paste0(indir,samples_list[i],path_postfix_filtered_matrix)
  obj.data <- Read10X(data.dir = datadir)
  this.obj <- CreateSeuratObject(counts = obj.data,
                                 # Don't accept data with fewer than 3 cells
                                 min.cells=3,
                                 # Don't accept data with fewer than 200 genes
                                 min.features=200)
  cells_pre_filter <- ncol(this.obj)
  cells_per_sample[nrow(cells_per_sample)+1,] <- c(obj.name,"pre_filter",
                                                   cells_pre_filter,0,0,0)
  
  #calculate the percentage of mitochondrial genes in each cell
  this.obj[["percent.mt"]] <- PercentageFeatureSet(this.obj, 
                                                   pattern = "^MT-")
  
  #Visualize QC metrics
  pdf(paste0(outdir,"/01_qc/",obj.name,"_pre_qc_plot.pdf"))
  print(VlnPlot(this.obj, features = c("nFeature_RNA",
                                       "nCount_RNA",
                                       "percent.mt"),
                ncol = 3))
  print(FeatureScatter(this.obj, feature1 = "nCount_RNA",
                       feature2 = "percent.mt"))
  print(FeatureScatter(this.obj, feature1 = "nCount_RNA",
                       feature2 = "nFeature_RNA"))
  dev.off()
  
  #filter based on QC cutoffs
  percent.mt.cutoff <- 5
  nfeature.min.cutoff <- 750
  nfeature.max.cutoff <- quantile(this.obj$nFeature_RNA, 0.99)
  this.obj <- subset(this.obj,subset = nFeature_RNA > nfeature.min.cutoff & 
                       nFeature_RNA < nfeature.max.cutoff & 
                       percent.mt < percent.mt.cutoff )
  
  # Visualize post QC metrics
  pdf(paste0(outdir,"/01_qc/",obj.name,"_post_qc_plot.pdf"))
  print(VlnPlot(this.obj, features = c("nFeature_RNA",
                                       "nCount_RNA",
                                       "percent.mt"),
                ncol = 3))
  print(FeatureScatter(this.obj, feature1 = "nCount_RNA",
                       feature2 = "percent.mt"))
  print(FeatureScatter(this.obj, feature1 = "nCount_RNA",
                       feature2 = "nFeature_RNA"))
  dev.off()
  
  #add metadata
  this.obj$sample_name <- substr(obj.name,1,9)
  this.obj$library_name <- obj.name
  this.obj$data_type <- "single-cell"
  this.obj$species <- "human"
  this.obj$prv_used <- "no"
  this.obj$library_prep_batch <- "batch_2"
  this.obj$experiment <- "new"
  
  print(paste0("Metadata for sample ",samples_list[i]," post-QC:"))
  print(head(this.obj))
  
  #get number of cells post filtering
  cells_post_filter <- ncol(this.obj)
  cells_per_sample[nrow(cells_per_sample)+1,] <- c(obj.name,
                                                   "post_filter",
                                                   cells_post_filter,
                                                   nfeature.min.cutoff,
                                                   nfeature.max.cutoff,
                                                   percent.mt.cutoff)
  
  #print the results to a text file
  sink(file.path(outdir,"01_qc/UL_transcripts_post_qc_sample_8_grch38.txt"), append = TRUE)
  genes_of_interest <- c("UL15","UL18","UL19","UL20","UL21","UL25","mRFP1","YFP","ChR2(H134R)")
  print(paste0("Sample ",samples_list[i], " :"))
  print(rowSums(this.obj[genes_of_interest,]@assays$RNA@data))
  sink()
  
  #save the filtered Seurat object
  saveRDS(this.obj, file = paste0(outdir,"/01_qc/",obj.name,".rds"))
  
}

#save the qc metadata
write.csv(cells_per_sample, 
          file = paste0(outdir,"/01_qc/qc_per_sample_metadata_sample_8_grch38.csv"), 
          row.names = FALSE)

#create a plot for cell counts before and after filtering for each sample
pdf(paste0(outdir, "/01_qc/","qc_cells_per_sample_plot_sample_8_grch38.pdf"),
    width = 10)
print(ggplot(cells_per_sample, 
             aes(factor(sample), 
                 as.numeric(number_of_cells), 
                 fill = rev(filter_stage))) + 
        geom_bar(stat="identity", 
                 position = "dodge") + 
        scale_fill_brewer(palette = "Set1",
                          labels = c("Pre-QC", "Post-QC"))+
        guides(fill=guide_legend(title="QC stage")) +
        xlab("Sample id") +
        ylab("Number of cells") +
        ggtitle("Cells per sample - Pre and Post QC") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
) 
dev.off()

#save the session info
writeLines(capture.output(sessionInfo()), 
           file.path(outdir, "sessionInfo_step01_qc_sample_8_grch38.txt"))

#record logs
print("***** Script completed! *****")

########################## END ##########################
