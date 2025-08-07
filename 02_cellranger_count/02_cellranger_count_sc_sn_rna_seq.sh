#!/bin/bash
#$ -cwd
#$ -o /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/tmp/
#$ -e /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/tmp/
#$ -pe smp 4
#$ -l mem_free=50G
#$ -l scratch=50G
#$ -l h_rt=10:00:00
#$ -j yes

###############################################################################
## Project ID: GB-LZ-1373
## Authors: Ayushi Agrawal
##
## Script Goal: run Cell Ranger count
##
## Expected command line arguments:
##    args[1]: Name of the directory containing all outputs
##    args[2]: Path of the directory contianing the fastq files
##    args[3]: Sample name to be processed by Cell Ranger count
##    args[4]: Path to the Cell Ranger-compatible transcriptome reference
##    args[5]: Set to true to include intronic reads in count 
##    args[6]: Path of the output directory
##
## Example running this script on wynton:
##    qsub cellranger_count_sc_sn_rna_seq.sh \
##	  Sample1A data/fastq_files Sample1A assets/refdata-lz-1373-mRatBN7-chr2-yfp-prv \
##	  True results/01_cellranger_count
###############################################################################

#setup paths
base_dir=/gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023
container_dir=$base_dir/assets
export SINGULARITY_BINDPATH="$base_dir"

#assign the argument values to variables
sample_out_id=$1
fastq_dir=$2
sample_name=$3
transcriptome_dir=$4
include_intron_val=$5
outdir=$6

#chnage working directory to output directory
cd $outdir

#run cellranger count
if [ "$include_intron_val" == "False" ]; then
	singularity exec $container_dir/cellranger-v6.1.1.sif cellranger count \
          --id=$sample_out_id \
          --transcriptome=$transcriptome_dir \
          --fastqs=$fastq_dir \
          --sample=$sample_name
else
	singularity exec $container_dir/cellranger-v6.1.1.sif cellranger count \
	  --id=$sample_out_id \
	  --include-introns=$include_intron_val \
	  --transcriptome=$transcriptome_dir \
	  --fastqs=$fastq_dir \
	  --sample=$sample_name
fi


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"


################### END ###################