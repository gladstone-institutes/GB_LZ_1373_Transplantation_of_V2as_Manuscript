#!/bin/bash
#$ -cwd
#$ -o /gladstone/bioinformatics/projects/lz-1266-lana-zholudeva-deepak-srivastava-scrnaseq-hg38-jan-2022/tmp/
#$ -e /gladstone/bioinformatics/projects/lz-1266-lana-zholudeva-deepak-srivastava-scrnaseq-hg38-jan-2022/tmp/
#$ -pe smp 4
#$ -l mem_free=50G
#$ -l scratch=50G
#$ -l h_rt=10:00:00
#$ -j yes

base_dir=/gladstone/bioinformatics/projects/lz-1266-lana-zholudeva-deepak-srivastava-scrnaseq-hg38-jan-2022
fastq_dir=$base_dir/data/gi-LZ3672
container_dir=$base_dir/assets/containers
export SINGULARITY_BINDPATH="$base_dir"

#assign the argument values to variables
sample_name=$1
sample_number=$(awk -F'_' '{print $2}' <<< $sample_name)
transcriptome_dir=$2
outdir=$3

#chnage working directory to output directory
cd $outdir

#if the sample is single-nucleus sample, use include_introns argument
if [[ "$sample_number" =~ ^(01|02|03|04)$ ]]; then 
  singularity exec $container_dir/cellranger-v6.1.1.sif cellranger count \
  --id=$sample_name \
  --transcriptome=$transcriptome_dir \
  --fastqs=$fastq_dir \
  --sample=$sample_name
else 
  singularity exec $container_dir/cellranger-v6.1.1.sif cellranger count \
  --id=$sample_name \
  --include-introns=true \
  --transcriptome=$transcriptome_dir \
  --fastqs=$fastq_dir \
  --sample=$sample_name
fi

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"


################### END ###################
