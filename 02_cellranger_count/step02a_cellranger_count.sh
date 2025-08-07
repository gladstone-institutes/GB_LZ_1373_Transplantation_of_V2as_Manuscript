#!/bin/bash

basedir=/gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023

# path to the cell ranger script
cellranger_script=$basedir/scripts/GB_LZ_1373_Transplantation_of_V2as_Manuscript/02_cellranger_count/02_cellranger_count_sc_sn_rna_seq.sh

#make the results output directory
outdir=$basedir/results/02_cellranger_count
mkdir -p $outdir

#transcriptome to use
transcriptome_dir=$basedir/results/01_custom_reference_genome/reference_genomes/refdata-lz-1373-grch38-chr2-yfp-prv

# run cell rnager count for all samples (run #1)
# qsub cell ranger count script for all 5 samples 
for d in /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/data/gi-LZ3672/LZ-3672_0[1-4]*_R1*; do
	#get sample name from the file name
	smp_filename=$(basename $d)
	smp=$(awk -F'_S' '{print $1}' <<< $smp_filename)

	#parameters for cellranger count
	out_id=${smp}_human
	include_intron_val=False
	fastq_dir=$basedir/data/gi-LZ3672
	
	# submit job
	qsub $cellranger_script $out_id $fastq_dir $smp $transcriptome_dir $include_intron_val $outdir
done

# run cell ranger count script for all samples (run #2)
for d in /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/data/gi_LZ4087/LZ408{2,7}_0[3-8]*_L001_R1*; do 
	#get sample name from the file name
	smp_filename=$(basename $d)
	smp=$(awk -F'_S' '{print $1}' <<< $smp_filename)
	out_id=${smp}_human

	#get sample number from the sample name
	sample_number=$(awk -F'_' '{print $2}' <<< $smp)

	#path to the fastq files
	fastq_dir=$basedir/data/gi_LZ4087

	#based on the sample number, assign the include_introns values
	if [[ "$sample_number" =~ ^(08) ]]; then 
		include_intron_val=False
	else
		include_intron_val=True
	fi

	#submit job
	qsub $cellranger_script $out_id $fastq_dir $smp $transcriptome_dir $include_intron_val $outdir
done


################### END ###################