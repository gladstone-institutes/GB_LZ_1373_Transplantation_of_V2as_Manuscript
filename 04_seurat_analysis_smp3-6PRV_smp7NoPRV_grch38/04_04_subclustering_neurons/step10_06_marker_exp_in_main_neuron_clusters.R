#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-LZ-1373
## Authors: Ayushi Agrawal
## Script Goal: Investigate three approaches to identify neurons for 
##              sub clustering: 
##                (1) RBFOX3 > 0 expression
##                (2) exact cluster IDs
##                (3) GFAP ≤ 0 expression.
###############################################################################

# Setup -----------------------------------------------------------
# Load required packages
library(dplyr)
library(Seurat)
library(ggplot2)

# define variables
basedir <- "~/Dropbox (Gladstone)/GB-LZ-1373/results"
input <- file.path(basedir, 
                   paste0("09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/",
                          "09_clustering_cell_type_no_CB_no_DF_HumanCells/",
                          "20PC_0.6res_rmBatchEffect_rmCluster16/",
                          "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_",
                          "rmBatchEffect_rmCluster16_clustered_and_cell_typed.rds"))
output_dir <- file.path(basedir,
                        "10_subcluster_neurons_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF",
                        "03_marker_exp_in_main_neuron_clusters")
selected_neuron_clusters <- c(1, 6, 10, 12, 15, 16)

# create the results folders
if (!(dir.exists(output_dir))) {
  dir.create(output_dir,recursive = TRUE)
}

# set working directory
setwd(output_dir)

#set seed
set.seed(42)

#increase RAM usage so the futurized Seurat functions... 
#can access larger global variables
options(future.globals.maxSize = 4000 * 1024^2)

# save outputs to a file
sink("marker_exp_in_main_neuron_clusters.txt")



# Main analysis -----------------------------------------------------------
# read in the data
data <- readRDS(input)

# check distrbution of clusters
print("Main clusters: ")
print(table(data$seurat_clusters))

# Question: Can we use RBFOX3 expression to identify neurons from the whole data?
# check RBFOX3 expression in the whole data
RBFOX3_expr <- FetchData(object = data, vars = "RBFOX3")
RBFOX3_expr_data <- data[, which(x = rowSums(RBFOX3_expr) > 0)]
print("RBFOX3>0 expressing cells across all clusters:")
print(table(RBFOX3_expr_data$seurat_clusters))
print(paste0("**Note: There are some cells in almost all clusters that express RBFOX3."))
print(paste0("**Lana's reply: RFBOX3 that is expressed in other than neuronal clusters",
             "are at the level of noise especially when checking for other neuronal ",
             "genes against the entire dataset. So, filter down to keep cells with both ",
             "RBFOX3>0 and those that are in clusters 1,6,10,12,15 and 16 for neuron sub-clsutering."))

pdf("FeaturePlot_neuron_and_glial_all_data_RBFOX3.pdf")
print(FeaturePlot(data, "RBFOX3", order = T, raster = F))
dev.off()

# Subset the data to include only main neuron clusters
print(paste0("Neuron clusters identified by Lana: 1, 6, 10, 12, 15, and 16"))
data_neuron_subset <- subset(data, idents = selected_neuron_clusters)

# filter cell will RBFOX3 
expr <- FetchData(object = data_neuron_subset, vars = "RBFOX3")
data_neuron_subset_rbfox3 <- data_neuron_subset[, which(x = rowSums(expr) > 0)]
percent_cells_rbfox3 <- round((ncol(data_neuron_subset_rbfox3)/ncol(data_neuron_subset))*100, 2)
print(paste0("There are a total of ", ncol(data_neuron_subset), 
             " nuclei in neuron clusters 1, 6, 10, 12, 15, and 16 combined. ",
             "When we filter these nuclei for RBFOX3 normalized expression > 0, ",
             "only ", ncol(data_neuron_subset_rbfox3), " (i.e. ", percent_cells_rbfox3,
             "%) nuclei remain for sub-clsutering."))

# check the distribution of RBFOX3 expression in data_neuron_subset between
# nuclei expressing GFAP and those not expressing GFAP
# Check if "GFAP" is in the active assay and normalized counts are present
# Add a metadata column for GFAP expression
data_neuron_subset$GFAP_meta <- ifelse(FetchData(data_neuron_subset, vars = "GFAP") > 0, 1, 0)

# Extract data for plotting
plot_data <- FetchData(data_neuron_subset, vars = c("RBFOX3", "GFAP_meta"))

# Convert GFAP_meta to a factor for better plotting
plot_data$GFAP_meta <- factor(plot_data$GFAP_meta, labels = c("GFAP <= 0", "GFAP > 0"))

# Create a violin plot
pdf("RBFOX3_expression_distribution_in_GFAP_expressing_neurons.pdf")
print(ggplot(plot_data, aes(x = GFAP_meta, y = RBFOX3, fill = GFAP_meta)) +
  geom_violin(trim = FALSE, scale = "width") + # Violin plot
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") + # Add boxplot for summary
  labs(title = "RBFOX3 expression in neuron clusters 1, 6, 10, 12, 15, and 16",
       x = "GFAP expression in cells",
       y = "RBFOX3 Expression") +
  theme_bw() +
  scale_fill_manual(values = c("skyblue", "lightpink")))
dev.off()
print(paste0("Looked at the violin plot distribution for RBFOX3 expression in neuron ",
             "clusters 1, 6, 10, 12, 15, and 16. RBFOX3 expression is zero inflated with",
             " many cells having 0 expression of this gene. Additionally, the distribution ",
             "of RBFOX3 in GFAP expressing nuclei is the same as that in nuclei not ",
             "expressing GFAP."))


# close the output file
sink()

# Write session info
writeLines(
  capture.output(sessionInfo()),
  "sessionInfo.txt"
)

print("***** Script Complete *****")



# END -----------------------------------------------------------