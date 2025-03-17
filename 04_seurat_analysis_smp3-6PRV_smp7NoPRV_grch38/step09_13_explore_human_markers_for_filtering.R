#!/usr/bin/env Rscript
###############################################################################
## Project ID: GB-LZ-1373
## Authors: Ayushi Agrawal
## Script Goal: Number of nuclei with different human cell markers as filters
## This script was run locally on Ayushi's laptop
###############################################################################

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)

#set all directory paths
basedir <- "~/Dropbox (Gladstone)/GB-LZ-1373/results/"
indir <- paste0(basedir,"/01_cellranger_count/")
path_postfix_filtered_matrix <- "/filtered_feature_bc_matrix"


#get the list of all rat/human spinal cord samples
samples_list <- c("LZ4082_03A_human","LZ4082_03B_human",
                  "LZ4082_04A_human","LZ4082_04B_human",
                  "LZ4087_05A_human","LZ4087_05B_human",
                  "LZ4087_06A_human","LZ4087_06B_human",
                  "LZ4087_07A_human_noPRV","LZ4087_07B_human_noPRV")

#create a data frame to record cell count before and after filtering per sample
cells_per_sample <- data.frame(sample=character(), 
                               filter_stage=character(), 
                               number_of_cells=numeric(),
                               nFeature_RNA_min_cutoff=numeric(),
                               nFeature_RNA_max_cutoff=numeric(),
                               percent.mt_cutoff=numeric(), 
                               stringsAsFactors = FALSE)

countsdf <- data.frame()

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
  cells_pre_QC <- ncol(this.obj)
  cells_per_sample[nrow(cells_per_sample)+1,] <- c(obj.name,"pre_QC",
                                                   cells_pre_QC,0,0,0)
  
  #calculate the percentage of mitochondrial genes in each cell
  this.obj[["percent.mt"]] <- PercentageFeatureSet(this.obj, 
                                                   pattern = "^MT-")
  
  #filter based on QC cutoffs
  percent.mt.cutoff <- 0.20
  nfeature.min.cutoff <- 200
  nfeature.max.cutoff <- quantile(this.obj$nFeature_RNA, 0.99)
  this.obj <- subset(this.obj,subset = nFeature_RNA > nfeature.min.cutoff & 
                       nFeature_RNA < nfeature.max.cutoff & 
                       percent.mt < percent.mt.cutoff )
  
  # get number of cells post filtering
  cells_post_QC <- ncol(this.obj)
  cells_per_sample[nrow(cells_per_sample)+1,] <- c(obj.name,
                                                   "post_QC",
                                                   cells_post_QC,
                                                   nfeature.min.cutoff,
                                                   nfeature.max.cutoff,
                                                   percent.mt.cutoff)
  
  
  # filter for human cells using Y chr genes and marker genes
  expr <- FetchData(object = this.obj, vars = c("USP9Y","UTY","RPS4Y1","YFP","ChR2(H134R)", "XRCC5"), slot = "counts")
  this.obj <- this.obj[, which(x = rowSums(expr) > 0)]
  expr <- FetchData(object = this.obj, vars = c("USP9Y","UTY","RPS4Y1","YFP","ChR2(H134R)", "XRCC5"), slot = "counts")
  print(as.data.frame(lapply(expr, function(c)sum(c!=0))))
  if(i >8 ){
    x <- as.data.frame(lapply(expr, function(c)sum(c!=0)))
    x$YFP <- 0
    x <- x[,c(1:3,6,4:5)]
    countsdf <- structure(rbind(countsdf,x), .Names = colnames(x))
  }else{
    countsdf <- structure(rbind(countsdf,as.data.frame(lapply(expr, function(c)sum(c!=0)))), .Names = colnames(as.data.frame(lapply(expr, function(c)sum(c!=0)))))
  }
  
  #get number of cells post filtering
  cells_post_gene_filter <- ncol(this.obj)
  cells_per_sample[nrow(cells_per_sample)+1,] <- c(obj.name,"post_YChrGene_filter",
                                                   cells_post_gene_filter,0,0,0)
}

# get the total number of nuclei
sum(as.numeric(cells_per_sample$number_of_cells[cells_per_sample$filter_stage == "post_YChrGene_filter"]))

# get the number of nuclei per marker gene
as.data.frame(lapply(countsdf, sum))


########################## END ##########################