#!/bin/bash
#$ -cwd
#$ -o /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/tmp/
#$ -e /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/tmp/
#$ -pe smp 1
#$ -l mem_free=40G
#$ -l scratch=50G
#$ -l h_rt=06:00:00
#$ -j yes

#file path variables
base_dir=/gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023
script_dir=$base_dir/scripts/GB-LZ-1373/08_seurat_analysis_samples_8_and_oldsample_grch38_no_cellbender_no_doubletfinder
container_dir=$base_dir/assets/containers
export APPTAINER_BINDPATH="$base_dir"

apptainer exec $container_dir/r_seurat_4.1_harmony_gb_lz_1266.sif Rscript $script_dir/08_01_seurat_qc_sample_8_grch38_no_cellbender_no_doubletfinder.R
  
## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
  
################### END ###################
