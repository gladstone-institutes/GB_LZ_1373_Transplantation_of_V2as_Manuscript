#!/bin/bash
#$ -cwd
#$ -o /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/tmp/
#$ -e /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/tmp/
#$ -pe smp 1
#$ -l mem_free=60G
#$ -l scratch=50G
#$ -l h_rt=04:00:00
#$ -j yes

#file path variables
base_dir=/gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023
script_dir=$base_dir/scripts/GB-LZ-1373/10_subcluster_neurons_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF
container_dir=$base_dir/assets/containers
export APPTAINER_BINDPATH="$base_dir"

apptainer exec $container_dir/seurat-4-3_soupx-1-6-1_doubletfinder_67fb8b5.sif Rscript \
$script_dir/10_01_subset_visualize.R \
--input_processed $base_dir/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/15_clustering_cell_type_no_CB_no_DF_HumanCells_NoPRV/20PC_0.6res_rmBatchEffect_rmCluster15/smp3-7_grch38_noCB_noDF_HumanCells_rmBatchEffect_NoPRV_20PC_0.6res_rmCluster15_clustered_and_cell_typed.rds \
--input_orig $base_dir/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/14_merge_and_visualize_no_CB_no_DF_HumanCells_NoPRV/smp3-7_grch38_noCB_noDF_HumanCells_NoPRV_merged_data.rds \
--output_dir $base_dir/results/10_subcluster_neurons_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/04_subset_and_visualize_NoPRV \
--output_prefix "neurons_smp3-7_grch38_noCB_noDF_NoPRV" \
--npcs 30 \
--cluster_subset "2,5,10,11,13,16,17,18"


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

################### END ###################
