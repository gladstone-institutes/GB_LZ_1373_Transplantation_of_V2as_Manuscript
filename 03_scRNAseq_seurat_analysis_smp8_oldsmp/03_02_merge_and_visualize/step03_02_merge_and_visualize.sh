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
data_dir=/gladstone/bioinformatics/projects
base_dir=$data_dir/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023
script_dir=$base_dir/scripts/GB_LZ_1373_Transplantation_of_V2as_Manuscript/03_scRNAseq_seurat_analysis_smp8_oldsmp
container_dir=$base_dir/assets
export APPTAINER_BINDPATH="$data_dir"

apptainer exec $container_dir/r_seurat_4.1_harmony_gb_lz_1266.sif Rscript \
$script_dir/03_02_merge_and_visualize/03_02_merge_and_visualize.R
  
## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
  
################### END ###################