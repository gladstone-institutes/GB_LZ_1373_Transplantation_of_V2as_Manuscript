#!/usr/bin/env Rscript
#run locally on Ayushi's laptop

#load required packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)

#set all directory paths
basedir <- paste0("~/Dropbox (Gladstone)/GB-LZ-1373/results/",
                  "10_subcluster_neurons_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/02_neuron_subclustering/")
setwd(basedir)

dir_list <- list.dirs(".", full.names = F, recursive = F)

for(i in 1:length(dir_list)){
  #load the Seurat object --------------------------------------------------------
  file_name <- list.files(dir_list[i], ".rds")
  data <- readRDS(file.path(dir_list[i], file_name))
  
  
  # generate heatmap -------------------------------------------------------------
  # heatmap for top 10 marker genes for each cluster
  de_markers <- read.csv(file.path(dir_list[i], 
                                   list.files(dir_list[i], 
                                              "_FindAllMarkers_marker_genes_per_cluster.csv")))
  de_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
  pdf(file = file.path(dir_list[i], paste0("neuron_subclusters_", dir_list[i], "_heatmap.pdf")),
      width = 28,
      height = 18)
  print(DoHeatmap(data, features = top10$gene) + 
          scale_fill_gradientn(colors = c("#33BBEE", "white", "#EE3377")) + 
          ggtitle(paste0("Heatmap of marker genes for neuron subclustering at ",dir_list[i], 
                         " \n(single-nucleus RNA-seq human cells)")))
  dev.off()
  
  
  # dotplots --------------------------------------------------------------------- 
  features_list <- list(general_population = c("SOX2", "PAX6", "MKI67", "slc1a3", "pax6",
                                                   "GLIS3", "CSPG4", "PLP1", "PDGFRA", "Olig3",
                                                   "MOG", "MAG", "MPZ", "PMP22"),
                            astrocytes = c("GFAP", "S100B", "AQP4", "slc1a3", "pax6", "GABBR2",
                                           "SLC6A11", "ADAMTS9", "PDE10A", "CPAMD8", "STAMBP",
                                           "Slc1a2", "Slc1a3", "Slc6a1", "Slc6a11", "Slc7a10",
                                           "Kcnj10", "Kcnj16"),
                            neurons = c("TUBB3", "MAP2", "NRG1", "GRM1", "SLC17A6", 
                                        "ROBO1", "ROBO2", "BCAS1", "TAC1", "SOX5", "PDE11A",
                                        "GAD1", "GAD2", "PAX2", "SLC6A5", "ADARB2", "PDYN", 
                                        "RORB", "ZFHX3", "NFIX", "MAF", "MAFA", "PAM", "TAC1",
                                        "TAC3", "NMU"))
  # None of the requested variables were found: TBFOX3
  # generate the dot plots
  for(f in 1:length(features_list)){
    features_to_plot <- unique(toupper(features_list[[f]]))
    pdf(file = file.path(dir_list[i], paste0("neuron_subclusters_", dir_list[i], "_",
                                             names(features_list)[f], "_dotplot.pdf")),
        width = 10)
    print(DotPlot(data, features = features_to_plot, 
                  cols = c("#33BBEE", "#EE3377")) + 
            RotatedAxis() + coord_flip())
    dev.off()
  }
  
  
  # >0 expression for percentage calculation for all marker genes ----------------
  percent_exp <- data.frame(gene_id = character(),
                            percentage_of_cells_with_greater_than_zero_expression = numeric(),
                            stringsAsFactors = F)
  for(f in unique(toupper(unlist(features_list)))){
    x <- FetchData(data, vars = f)[,1]
    high_exp_cells_percent <- round(((length(x[x>0])*100)/length(x)),2)
    percent_exp[nrow(percent_exp)+1,] <- c(f,high_exp_cells_percent)
  }
  write.csv(percent_exp,
            file.path(dir_list[i], paste0("neuron_subclusters_", dir_list[i], 
                                          "_percentage_expression_of_marker_genes.csv")),
            row.names = F)
  
  
  # cells per sample per cluster -------------------------------------------------
  ##the number of cells in each mouse in each cluster
  all_metadata <- data[[]]
  
  #make a table of the required data
  smp_cluster_counts <- unique(data@meta.data %>%
                                 group_by(sample_name) %>%
                                 mutate(total_numbers_of_cells_per_sample = n()) %>%
                                 group_by(seurat_clusters, .add=TRUE) %>%
                                 mutate(number_of_cells_per_sample_in_cluster = n(), 
                                        proportion_of_cells_per_sample_in_cluster =
                                          number_of_cells_per_sample_in_cluster/total_numbers_of_cells_per_sample,
                                        percent_of_cells_per_sample_in_cluster = 
                                          proportion_of_cells_per_sample_in_cluster*100 ) %>%
                                 select(sample_name,
                                        library_prep_batch,
                                        seurat_clusters,
                                        total_numbers_of_cells_per_sample,
                                        number_of_cells_per_sample_in_cluster, 
                                        proportion_of_cells_per_sample_in_cluster,
                                        percent_of_cells_per_sample_in_cluster))
  colnames(smp_cluster_counts)[1:3] <- c("SampleName","LibraryPrepBatch","SeuratClusterID")
  smp_cluster_counts <- smp_cluster_counts[order(smp_cluster_counts$SampleName),]
  write.csv(smp_cluster_counts, 
            file = file.path(dir_list[i], paste0("neuron_subclusters_", dir_list[i], 
                                                 "_counts_per_sample_per_cluster.csv")),
            row.names = FALSE)
  
  
  
  # cells expressing UL18, UL19 and UL25 for each animal -------------------------
  ##the number of cells in each mouse in each cluster
  all_metadata <- data[[]]
  
  expr_df <- FetchData(data, vars = c("UL18", "UL19", "UL25"))
  
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
            file = file.path(dir_list[i], paste0("neuron_subclusters_", dir_list[i], 
                                                 "_cell_count_per_sample_with_UL_gene_exp_greater_than_zero.csv")),
            row.names = FALSE)
  
  
}


# END --------------------------------------------------------------------------