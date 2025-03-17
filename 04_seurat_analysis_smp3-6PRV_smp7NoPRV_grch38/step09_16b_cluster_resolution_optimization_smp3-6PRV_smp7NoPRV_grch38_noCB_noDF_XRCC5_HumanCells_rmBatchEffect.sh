#!/usr/bin/env bash
#$ -cwd
#$ -o /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/tmp/
#$ -e /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/tmp/
#$ -pe smp 8
#$ -l mem_free=40G
#$ -l scratch=50G
#$ -l h_rt=06:00:00
#$ -j yes

data_dir=/gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023
script_dir=$data_dir/scripts/GB-LZ-1373/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF
results_dir=$data_dir/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF
container_dir=$data_dir/assets/containers
export APPTAINER_BINDPATH="$data_dir"

singularity exec $container_dir/seurat-4-3_soupx-1-6-1_doubletfinder_67fb8b5.sif \
Rscript $script_dir/09_03_cluster_resolution_optimization.R \
    --input $results_dir/11_merge_and_visualize_noCB_noDF_XRCC5_HumanCells/smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_XRCC5_HumanCells_merged_data_sct_pca_umap.rds \
    --output $results_dir/12_cluster_resolution_optimization_noCB_noDF_XRCC5_HumanCells/15_PCs_rmBatchEffect \
    --output_prefix "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_XRCC5_HumanCells_rmBatchEffect_15PC" \
    --mode "scRNA" \
    --cores ${NSLOTS} \
    --ndim 15 \
    --metadata "sample_name" \
    --within_batch_var "library_prep_batch"

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

######################## END ########################