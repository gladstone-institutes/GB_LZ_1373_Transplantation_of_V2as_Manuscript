#!/usr/bin/env Rscript
###############################################################################
## Project ID: GB-LZ-1373
## Authors: Ayushi Agrawal
## Script Goal: Try label transfer for single-nucleus data. Previous experiment was
##              was used as reference data.
## Note: This script was run locally on Ayushi's laptop
## Example running this script on wynton:
##    Rscript step09_18a_transfer_label_cell_annotations.R
###############################################################################


# load required packages -------------------------------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)


# load data and set all paths --------------------------------------------------
reference <- readRDS(paste0("~/Dropbox (Gladstone)/GB-LZ-1266/results_rerun_on_wynton/",
                            "03_seurat_analysis_newQCcutoffs_grch38_chr2_yfp/data/",
                            "02_clustering/single_nucleus_data/20_PCs_0.6_resolution/",
                            "merged_data_newQCcutoffs_processed_single_nucleus_data_20_PCs_0.6_res.rds"))
query <- readRDS(paste0("~/Dropbox (Gladstone)/GB-LZ-1373/results/",
                        "09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/",
                        "09_clustering_cell_type_no_CB_no_DF_HumanCells/",
                        "20PC_0.6res_rmBatchEffect_rmCluster16/smp3-6PRV_smp7NoPRV",
                        "_grch38_noCB_noDF_HumanCells_rmBatchEffect_rmCluster16_clustered_and_cell_typed.rds"))
outdir <- paste0("~/Dropbox (Gladstone)/GB-LZ-1373/results/",
                 "09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/",
                 "09_clustering_cell_type_no_CB_no_DF_HumanCells/",
                 "20PC_0.6res_rmBatchEffect_rmCluster16")
setwd(outdir)


# set cell type labels for reference -------------------------------------------
# using annotations in Fig 2A (https://www.biorxiv.org/content/10.1101/2024.01.11.575264v1.full)
reference$celltype <- ifelse(reference$seurat_clusters %in% c(16,19,7,3, 12, 8,9,5), "Ex",
                             ifelse(reference$seurat_clusters %in% c(13,4,10), "Inh",
                                    ifelse(reference$seurat_clusters %in% c(0,1,2,11), "GMAst",
                                           ifelse(reference$seurat_clusters %in% c(6), "WMAst",
                                                  ifelse(reference$seurat_clusters %in% c(15), "AstPgen",
                                                         ifelse(reference$seurat_clusters %in% c(14,18), "Oligo","OligoPgen"))))))


# label transfer ---------------------------------------------------------------
dat.anchors <- FindTransferAnchors(reference = reference, 
                                   query = query, 
                                   normalization.method = "SCT",
                                   reference.assay = "SCT",
                                   query.assay = "SCT",
                                   reference.reduction = "pca",
                                   dims = 1:30)
predictions <- TransferData(anchorset = dat.anchors, 
                            refdata = reference$celltype, 
                            dims = 1:30)
dat.query <- AddMetaData(query, metadata = predictions)

result_df <- as.data.frame.matrix(table(dat.query$predicted.id, dat.query$seurat_clusters))
colnames(result_df) <- paste0("Cluster_",colnames(result_df))

# Write the DataFrame to the file
write.csv(result_df, 
          "label_transfer_annotations_from_old_experiment_single_nucleus_data.csv")

# Note to add at the beginning of the file
note <- paste0("Note: Seurat's transfer label was used to annotate the single-cell data using",
               " annotations in Fig 2A (https://www.biorxiv.org/content/10.1101/2024.01.11.575264v1.full). ",
               "The rows are the annotations. The columns are the cluster ids in the new data. ",
               "The values represent the number of single-cells in a cluster predicted for the respective annotation.")

# Read the existing CSV content
csv_content <- readLines("label_transfer_annotations_from_old_experiment_single_nucleus_data.csv")

# Combine the note and the original CSV content
csv_with_note <- c(note, csv_content)

# Write the modified content back to the CSV file
writeLines(csv_with_note, "label_transfer_annotations_from_old_experiment_single_nucleus_data.csv")



# END --------------------------------------------------------------------------