#!/usr/bin/env Rscript
###############################################################################
## Project ID: GB-LZ-1373
## Authors: Ayushi Agrawal
## Script Goal: Try label transfer for single-cell data. Previous experiment was
##              was used as reference data.
## Note: This script was run locally on Ayushi's laptop
## Example running this script on wynton:
##    Rscript step08_07b_transfer_label_cell_annotations.R
###############################################################################


# load required packages -------------------------------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)


# load data and set all paths --------------------------------------------------
reference <- readRDS(paste0("~/Dropbox (Gladstone)/GB-LZ-1266/results_rerun_on_wynton/",
                            "03_seurat_analysis_newQCcutoffs_grch38_chr2_yfp/data/",
                            "02_clustering/single_cell_data/20_PCs_0.09_resolution/",
                            "merged_data_newQCcutoffs_processed_single_cell_data_20_PCs_0.09_res.rds"))
query <- readRDS(paste0("~/Dropbox (Gladstone)/GB-LZ-1373/results/",
                        "08_seurat_analysis_samples_8_and_oldsample_grch38_",
                        "no_cellbender_no_doubletfinder/04_clustering_cell_type/",
                        "15PC_0.1res/samples_8_and_oldsample_grch38_clustered_",
                        "and_cell_typed.rds"))
outdir <- paste0("~/Dropbox (Gladstone)/GB-LZ-1373/results/",
                 "08_seurat_analysis_samples_8_and_oldsample_grch38_no_cellbender_no_doubletfinder",
                 "/04_clustering_cell_type/15PC_0.1res")
setwd(outdir)


# set cell type labels for reference -------------------------------------------
# using annotations in Fig 2A (https://www.biorxiv.org/content/10.1101/2024.01.11.575264v1.full)
reference$celltype <- ifelse(reference$seurat_clusters == 0, "V2;p2",
                             ifelse(reference$seurat_clusters == 1, "Spinal Progenitors",
                                    ifelse(reference$seurat_clusters == 2, "V0;V1;V3",
                                           ifelse(reference$seurat_clusters == 3, "p3;glial/PNS progenitors","MNs"))))


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
          "label_transfer_annotations_from_old_experiment_single_cell_data.csv")

# Note to add at the beginning of the file
note <- paste0("Note: Seurat's transfer label was used to annotate the single-cell data using",
               " annotations in Fig 2A (https://www.biorxiv.org/content/10.1101/2024.01.11.575264v1.full). ",
               "The rows are the annotations. The columns are the cluster ids in the new data. ",
               "The values represent the number of single-cells in a cluster predicted for the respective annotation.")

# Read the existing CSV content
csv_content <- readLines("label_transfer_annotations_from_old_experiment_single_cell_data.csv")

# Combine the note and the original CSV content
csv_with_note <- c(note, csv_content)

# Write the modified content back to the CSV file
writeLines(csv_with_note, "label_transfer_annotations_from_old_experiment_single_cell_data.csv")



# END --------------------------------------------------------------------------