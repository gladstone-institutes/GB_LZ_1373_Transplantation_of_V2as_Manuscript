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

data_dir=/gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023
script_dir=$data_dir/scripts/GB_LZ_1373_Transplantation_of_V2as_Manuscript/03_scRNAseq_seurat_analysis_smp8_oldsmp
results_dir=$data_dir/results/03_scRNAseq_seurat_analysis_smp8_oldsmp
container_dir=$data_dir/assets
export APPTAINER_BINDPATH="$data_dir"

resolutions=(0.02 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2)
res=${resolutions[$SGE_TASK_ID-1]}

echo "Running task ID $SGE_TASK_ID with resolution = $res"


singularity exec $container_dir/seurat-4-3_soupx-1-6-1_doubletfinder_67fb8b5.sif \
Rscript $script_dir/03_03_clustering_cell_type/03_03_clustering_cell_type.R \
    --input $results_dir/02_merge_and_visualize/merged_data_samples_8_and_oldsample_grch38_post_sct_processed.rds \
    --output $results_dir/03_clustering_cell_type/15PC_${res}res \
    --output_prefix "samples_8_and_oldsample_grch38" \
    --tissue "Brain" \
    --ndim 15 \
    --resolution $res


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"


######################## END ########################