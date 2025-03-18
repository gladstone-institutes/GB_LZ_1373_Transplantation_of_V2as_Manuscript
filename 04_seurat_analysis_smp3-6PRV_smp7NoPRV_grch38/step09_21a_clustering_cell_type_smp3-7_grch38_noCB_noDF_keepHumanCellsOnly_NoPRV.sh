#!/usr/bin/env bash
#$ -cwd
#$ -o /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/tmp/
#$ -e /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/tmp/
#$ -pe smp 1
#$ -l mem_free=20G
#$ -l scratch=50G
#$ -l h_rt=02:00:00
#$ -j yes
#$ -t 1-16

# Define directories
data_dir=/gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023
script_dir=$data_dir/scripts/GB-LZ-1373/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF
results_dir=$data_dir/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF
container_dir=$data_dir/assets/containers
export APPTAINER_BINDPATH="$data_dir"

# List of resolutions to iterate over
resolutions=(0.02 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2)

# Get the resolution for this task
res=${resolutions[$SGE_TASK_ID-1]}

# Run clustering
singularity exec $container_dir/seurat-4-3_soupx-1-6-1_doubletfinder_67fb8b5.sif \
Rscript $script_dir/09_04_clustering_cell_type.R \
--input $results_dir/14_merge_and_visualize_no_CB_no_DF_HumanCells_NoPRV/smp3-7_grch38_noCB_noDF_HumanCells_NoPRV_merged_data_sct_pca_umap.rds \
--output $results_dir/15_clustering_cell_type_no_CB_no_DF_HumanCells_NoPRV/20PC_${res}res_rmBatchEffect \
--output_prefix "smp3-7_grch38_noCB_noDF_HumanCells_rmBatchEffect_NoPRV_20PC_${res}res" \
--tissue "Brain" \
--ndim 20 \
--resolution $res \
--batch_var "library_prep_batch"


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

######################## END ########################
