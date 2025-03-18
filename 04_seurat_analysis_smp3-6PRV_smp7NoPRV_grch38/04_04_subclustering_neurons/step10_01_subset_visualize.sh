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
--input $base_dir/results/09_seurat_analysis_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/09_clustering_cell_type_no_CB_no_DF_HumanCells/20PC_0.6res_rmBatchEffect_rmCluster16/smp3-6PRV_smp7NoPRV_grch38_noCB_noDF_HumanCells_rmBatchEffect_rmCluster16_clustered_and_cell_typed.rds  \
--output_dir $base_dir/results/10_subcluster_neurons_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF/01_subset_and_visualize \
--output_prefix "neurons_smp3-6PRV_smp7NoPRV_grch38_noCB_noDF" \
--npcs 20 \
--cluster_subset "1,6,10,12,15,16"

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

################### END ###################
