#!/usr/bin/env bash
#$ -cwd
#$ -o /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/tmp/
#$ -e /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/tmp/
#$ -pe smp 1
#$ -l mem_free=20G
#$ -l scratch=50G
#$ -l h_rt=02:00:00
#$ -j yes

data_dir=/gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023
script_dir=$data_dir/scripts/GB-LZ-1373/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF
results_dir=$data_dir/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF
container_dir=$data_dir/assets/containers
export APPTAINER_BINDPATH="$data_dir"


singularity exec $container_dir/seurat-4-3_soupx-1-6-1_doubletfinder_67fb8b5.sif \
Rscript $script_dir/09_04_clustering_cell_type.R \
    --input $results_dir/07_merge_and_visualize_no_CB_no_DF_HumanCells/smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_merged_data_sct_pca_umap.rds \
    --output $results_dir/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.6res_rmBatchEffect \
    --output_prefix "smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_rmBatchEffect" \
    --tissue "Brain" \
    --ndim 20 \
    --resolution 0.6 \
    --batch_var "library_prep_batch"


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

######################## END ########################