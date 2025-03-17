#!/usr/bin/env bash
#$ -cwd
#$ -o /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/tmp/
#$ -e /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/tmp/
#$ -pe smp 1
#$ -l mem_free=40G
#$ -l scratch=50G
#$ -l h_rt=20:00:00
#$ -j yes

data_dir=/gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023
script_dir=$data_dir/scripts/GB-LZ-1373/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF
results_dir=$data_dir/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF
container_dir=$data_dir/assets/containers
export APPTAINER_BINDPATH="$data_dir"


singularity exec $container_dir/seurat-4-3_soupx-1-6-1_doubletfinder_67fb8b5.sif \
Rscript $script_dir/09_04_clustering_cell_type.R \
    --input $results_dir/11_merge_and_visualize_noCB_noDF_XRCC5_HumanCells/smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_XRCC5_HumanCells_merged_data_sct_pca_umap.rds \
    --output $results_dir/13_clustering_cell_type_noCB_noDF_XRCC5_HumanCells/15PC_0.02res_rmBatchEffect \
    --output_prefix "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_XRCC5_HumanCells_rmBatchEffect" \
    --tissue "Brain" \
    --ndim 15 \
    --resolution 0.02 \
    --batch_var "library_prep_batch"

singularity exec $container_dir/seurat-4-3_soupx-1-6-1_doubletfinder_67fb8b5.sif \
Rscript $script_dir/09_04_clustering_cell_type.R \
    --input $results_dir/11_merge_and_visualize_noCB_noDF_XRCC5_HumanCells/smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_XRCC5_HumanCells_merged_data_sct_pca_umap.rds \
    --output $results_dir/13_clustering_cell_type_noCB_noDF_XRCC5_HumanCells/15PC_0.04res_rmBatchEffect \
    --output_prefix "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_XRCC5_HumanCells_rmBatchEffect" \
    --tissue "Brain" \
    --ndim 15 \
    --resolution 0.04 \
    --batch_var "library_prep_batch"

singularity exec $container_dir/seurat-4-3_soupx-1-6-1_doubletfinder_67fb8b5.sif \
Rscript $script_dir/09_04_clustering_cell_type.R \
    --input $results_dir/11_merge_and_visualize_noCB_noDF_XRCC5_HumanCells/smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_XRCC5_HumanCells_merged_data_sct_pca_umap.rds \
    --output $results_dir/13_clustering_cell_type_noCB_noDF_XRCC5_HumanCells/15PC_0.06res_rmBatchEffect \
    --output_prefix "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_XRCC5_HumanCells_rmBatchEffect" \
    --tissue "Brain" \
    --ndim 15 \
    --resolution 0.06 \
    --batch_var "library_prep_batch"

singularity exec $container_dir/seurat-4-3_soupx-1-6-1_doubletfinder_67fb8b5.sif \
Rscript $script_dir/09_04_clustering_cell_type.R \
    --input $results_dir/11_merge_and_visualize_noCB_noDF_XRCC5_HumanCells/smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_XRCC5_HumanCells_merged_data_sct_pca_umap.rds \
    --output $results_dir/13_clustering_cell_type_noCB_noDF_XRCC5_HumanCells/15PC_0.1res_rmBatchEffect \
    --output_prefix "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_XRCC5_HumanCells_rmBatchEffect" \
    --tissue "Brain" \
    --ndim 15 \
    --resolution 0.1 \
    --batch_var "library_prep_batch"

singularity exec $container_dir/seurat-4-3_soupx-1-6-1_doubletfinder_67fb8b5.sif \
Rscript $script_dir/09_04_clustering_cell_type.R \
    --input $results_dir/11_merge_and_visualize_noCB_noDF_XRCC5_HumanCells/smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_XRCC5_HumanCells_merged_data_sct_pca_umap.rds \
    --output $results_dir/13_clustering_cell_type_noCB_noDF_XRCC5_HumanCells/15PC_0.2res_rmBatchEffect \
    --output_prefix "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_XRCC5_HumanCells_rmBatchEffect" \
    --tissue "Brain" \
    --ndim 15 \
    --resolution 0.2 \
    --batch_var "library_prep_batch"

singularity exec $container_dir/seurat-4-3_soupx-1-6-1_doubletfinder_67fb8b5.sif \
Rscript $script_dir/09_04_clustering_cell_type.R \
    --input $results_dir/11_merge_and_visualize_noCB_noDF_XRCC5_HumanCells/smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_XRCC5_HumanCells_merged_data_sct_pca_umap.rds \
    --output $results_dir/13_clustering_cell_type_noCB_noDF_XRCC5_HumanCells/15PC_0.4res_rmBatchEffect \
    --output_prefix "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_XRCC5_HumanCells_rmBatchEffect" \
    --tissue "Brain" \
    --ndim 15 \
    --resolution 0.4 \
    --batch_var "library_prep_batch"

singularity exec $container_dir/seurat-4-3_soupx-1-6-1_doubletfinder_67fb8b5.sif \
Rscript $script_dir/09_04_clustering_cell_type.R \
    --input $results_dir/11_merge_and_visualize_noCB_noDF_XRCC5_HumanCells/smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_XRCC5_HumanCells_merged_data_sct_pca_umap.rds \
    --output $results_dir/13_clustering_cell_type_noCB_noDF_XRCC5_HumanCells/15PC_0.6res_rmBatchEffect \
    --output_prefix "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_XRCC5_HumanCells_rmBatchEffect" \
    --tissue "Brain" \
    --ndim 15 \
    --resolution 0.6 \
    --batch_var "library_prep_batch"

singularity exec $container_dir/seurat-4-3_soupx-1-6-1_doubletfinder_67fb8b5.sif \
Rscript $script_dir/09_04_clustering_cell_type.R \
    --input $results_dir/11_merge_and_visualize_noCB_noDF_XRCC5_HumanCells/smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_XRCC5_HumanCells_merged_data_sct_pca_umap.rds \
    --output $results_dir/13_clustering_cell_type_noCB_noDF_XRCC5_HumanCells/15PC_0.7res_rmBatchEffect \
    --output_prefix "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_XRCC5_HumanCells_rmBatchEffect" \
    --tissue "Brain" \
    --ndim 15 \
    --resolution 0.7 \
    --batch_var "library_prep_batch"



## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

######################## END ########################