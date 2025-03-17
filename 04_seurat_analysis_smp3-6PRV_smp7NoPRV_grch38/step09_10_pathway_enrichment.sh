#!/usr/bin/env bash

# Wrapper to run pseudobulk DE analysis on all cells
Rscript 09_10_pathway_enrichment.R \
--input '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.02res' \
--output '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.02res' \
--pathway_db '~/Dropbox (Gladstone)/GB-LZ-1373/assets/Hs_pathway_20240615.RData' \
--species 'human' \
--p_val_cutoff 0.01 \
--fold_change_cutoff 1


Rscript 09_10_pathway_enrichment.R \
--input '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.04res' \
--output '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.04res' \
--pathway_db '~/Dropbox (Gladstone)/GB-LZ-1373/assets/Hs_pathway_20240615.RData' \
--species 'human' \
--p_val_cutoff 0.01 \
--fold_change_cutoff 1


Rscript 09_10_pathway_enrichment.R \
--input '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.5res' \
--output '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.5res' \
--pathway_db '~/Dropbox (Gladstone)/GB-LZ-1373/assets/Hs_pathway_20240615.RData' \
--species 'human' \
--p_val_cutoff 0.01 \
--fold_change_cutoff 1


Rscript 09_10_pathway_enrichment.R \
--input '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.3res' \
--output '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.3res' \
--pathway_db '~/Dropbox (Gladstone)/GB-LZ-1373/assets/Hs_pathway_20240615.RData' \
--species 'human' \
--p_val_cutoff 0.01 \
--fold_change_cutoff 1


Rscript 09_10_pathway_enrichment.R \
--input '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.6res' \
--output '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.6res' \
--pathway_db '~/Dropbox (Gladstone)/GB-LZ-1373/assets/Hs_pathway_20240615.RData' \
--species 'human' \
--p_val_cutoff 0.01 \
--fold_change_cutoff 1


Rscript 09_10_pathway_enrichment.R \
--input '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.2res_rmBatchEffect' \
--output '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.2res_rmBatchEffect' \
--pathway_db '~/Dropbox (Gladstone)/GB-LZ-1373/assets/Hs_pathway_20240615.RData' \
--species 'human' \
--p_val_cutoff 0.01 \
--fold_change_cutoff 1


Rscript 09_10_pathway_enrichment.R \
--input '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.6res_rmBatchEffect' \
--output '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.6res_rmBatchEffect' \
--pathway_db '~/Dropbox (Gladstone)/GB-LZ-1373/assets/Hs_pathway_20240615.RData' \
--species 'human' \
--p_val_cutoff 0.01 \
--fold_change_cutoff 1


Rscript 09_10_pathway_enrichment.R \
--input '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.02res_rmBatchEffect' \
--output '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.02res_rmBatchEffect' \
--pathway_db '~/Dropbox (Gladstone)/GB-LZ-1373/assets/Hs_pathway_20240615.RData' \
--species 'human' \
--p_val_cutoff 0.01 \
--fold_change_cutoff 1


Rscript 09_10_pathway_enrichment.R \
--input '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.04res_rmBatchEffect' \
--output '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.04res_rmBatchEffect' \
--pathway_db '~/Dropbox (Gladstone)/GB-LZ-1373/assets/Hs_pathway_20240615.RData' \
--species 'human' \
--p_val_cutoff 0.01 \
--fold_change_cutoff 1


Rscript 09_10_pathway_enrichment.R \
--input '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.6res_rmBatchEffect_rmCluster16' \
--output '~/Dropbox (Gladstone)/GB-LZ-1373/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.6res_rmBatchEffect_rmCluster16' \
--pathway_db '~/Dropbox (Gladstone)/GB-LZ-1373/assets/Hs_pathway_20240615.RData' \
--species 'human' \
--p_val_cutoff 0.01 \
--fold_change_cutoff 1


# END -----------------------------------------------------