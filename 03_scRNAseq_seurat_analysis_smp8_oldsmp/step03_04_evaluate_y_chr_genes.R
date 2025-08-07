#!/usr/bin/env Rscript
#run locally on Ayushi's laptop

#load required packages
library(Seurat)
library(dplyr)
library(readxl)
library(ggplot2)

#set all directory paths
basedir <- "~/Dropbox (Gladstone)/GB-LZ-1373/results/03_scRNAseq_seurat_analysis_smp8_oldsmp/02_merge_and_visualize/"
setwd(basedir)

#load the Seurat object
merged_sc_data <- readRDS("merged_data_samples_8_and_oldsample_grch38_post_sct_processed.rds")

y_chr_genes <- read_xlsx("~/Dropbox (Gladstone)/github/GB_LZ_1373_Transplantation_of_V2as_Manuscript/reference/Y_Chromosome_Specific_Genes.xlsx")
y_chr_genes <- unique(y_chr_genes$`Gene Symbol`)
y_chr_genes <- y_chr_genes[y_chr_genes %in% rownames(merged_sc_data)]


# >0 expression of Y chr specific genes ----------------------------------------
percent_exp <- data.frame(gene_id = character(),
                          percentage_of_cells_with_greater_than_zero_expression = numeric(),
                          stringsAsFactors = F)
for(i in y_chr_genes){
  x <- FetchData(merged_sc_data, vars = i, layer = "counts")[,1]
  high_exp_cells_percent <- round(((length(x[x>0])*100)/length(x)),2)
  percent_exp[nrow(percent_exp)+1,] <- c(i,high_exp_cells_percent)
}
write.csv(percent_exp,
          "merged_data_samples_8_and_oldsample_grch38_percentage_SCTcounts_of_Ychr_genes.csv",
          row.names = F)


# violin plots of Y chr genes --------------------------------------------------
pdf(file = "merged_data_samples_8_and_oldsample_grch38_Ychr_genes_SCTcounts_vlnplot.pdf",
    height = 10)
print(VlnPlot(merged_sc_data, 
              layer = "counts",
              features = y_chr_genes))
dev.off()

pdf(file = "merged_data_samples_8_and_oldsample_grch38_Ychr_genes_SCTcounts_split_vlnplot.pdf",
    height = 10)
print(VlnPlot(merged_sc_data, 
              layer = "counts",
              features = y_chr_genes,
              split.by = "library_prep_batch",
              group.by = "library_prep_batch"))
dev.off()

pdf(file = "merged_data_samples_8_and_oldsample_grch38_Ychr_genes_SCTnormalizedcounts_vlnplot.pdf",
    height = 10)
print(VlnPlot(merged_sc_data,
              features = y_chr_genes))
dev.off()


#################### END ####################