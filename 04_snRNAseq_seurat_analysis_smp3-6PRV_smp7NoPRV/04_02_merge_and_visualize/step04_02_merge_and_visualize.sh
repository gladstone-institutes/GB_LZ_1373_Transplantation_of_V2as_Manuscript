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
script_dir=$base_dir/scripts/GB_LZ_1373_Transplantation_of_V2as_Manuscript/04_snRNAseq_seurat_analysis_smp3-6PRV_smp7NoPRV
container_dir=$base_dir/assets
export APPTAINER_BINDPATH="$base_dir"

apptainer exec $container_dir/seurat-4-3_soupx-1-6-1_doubletfinder_67fb8b5.sif Rscript \
$script_dir/04_02_merge_and_visualize/04_02_merge_and_visualize.R \
--input_dir $base_dir/results/04_snRNAseq_seurat_analysis_smp3-6PRV_smp7NoPRV/01_qc \
--output_dir $base_dir/results/04_snRNAseq_seurat_analysis_smp3-6PRV_smp7NoPRV/02_merge_and_visualize \
--output_prefix "smp3-7_grch38_noCB_noDF_HumanCells_NoPRV" \
--project "LZ_1373" \
--npcs 20
  
## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
  
################### END ###################