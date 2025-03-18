#!/bin/bash

basedir=/gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023

#path to the fastq files
fastq_dir=$basedir/data/gi_LZ4087

# path to the cell ranger script
cellranger_script=$basedir/scripts/GB-LZ-1373/01_cellranger_count/cellranger_count_sc_sn_rna_seq.sh

#make the results output directory
outdir=$basedir/results/01_cellranger_count
mkdir -p $outdir

#transcriptome to use  
transcriptome_rat=$basedir/assets/reference_genomes/refdata-lz-1266-mRatBN7-chr2-yfp-prv

# run cell ranger count script for all 8 samples 
for d in /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/data/gi_LZ4087/*_L001_R1*; do 
	#get sample name from the file name
	smp_filename=$(basename $d)
	smp=$(awk -F'_S' '{print $1}' <<< $smp_filename)

	#get sample number from the sample name
	sample_number=$(awk -F'_' '{print $2}' <<< $smp)

	#based on the sample number, assign the ouptut id, transcriptome and include_introns values
	if [[ "$sample_number" =~ ^(03|04|05|06|07) ]]; then 
		echo $sample_number
		
		out_id=${smp}_rat
		transcriptome_dir=$transcriptome_rat
		include_intron_val=true
		
		#run cellranger
		qsub $cellranger_script $out_id $fastq_dir $smp $transcriptome_dir $include_intron_val $outdir
        fi
done

################### END ###################
