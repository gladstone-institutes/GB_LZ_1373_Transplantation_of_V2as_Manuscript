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
script_dir=$data_dir/scripts/GB-LZ-1373/08_seurat_analysis_samples_8_and_oldsample_grch38_no_cellbender_no_doubletfinder
results_dir=$data_dir/results/08_seurat_analysis_samples_8_and_oldsample_grch38_no_cellbender_no_doubletfinder
container_dir=$data_dir/assets/containers
export APPTAINER_BINDPATH="$data_dir"


singularity exec $container_dir/seurat-4-3_soupx-1-6-1_doubletfinder_67fb8b5.sif \
Rscript $script_dir/08_04_clustering_cell_type.R \
    --input $results_dir/02_merge_and_visualize/merged_data_samples_8_and_oldsample_grch38_post_sct_processed.rds \
    --output $results_dir/04_clustering_cell_type/15PC_0.1res \
    --output_prefix "samples_8_and_oldsample_grch38" \
    --tissue "Brain" \
    --ndim 15 \
    --resolution 0.1


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

######################## END ########################
