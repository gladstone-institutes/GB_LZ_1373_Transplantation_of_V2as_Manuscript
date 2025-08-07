#!/bin/bash

# One snRNA-seq sample (sample 7), which had not been exposed to PRV and with low cell recovery (<1,000 cells), 
# unexpectedly showed expression of some PRV genes despite the absence of viral treatment. 
# This was likely due to spurious alignment or background noise. To avoid potential artifacts, 
# this sample was reprocessed using the PRV-excluded reference. No PRV expression was detected in the  
# scRNA-seq samples, which were of higher quality and showed no alignment to PRV sequences as expected; 
# therefore, the PRV-inclusive reference was retained for those samples to maintain consistency across the dataset.

basedir=/gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023

#path to the fastq files
fastq_dir=$basedir/data/gi_LZ4087

# path to the cell ranger script
cellranger_script=$basedir/scripts/GB_LZ_1373_Transplantation_of_V2as_Manuscript/02_cellranger_count/02_cellranger_count_sc_sn_rna_seq.sh

#make the results output directory
outdir=$basedir/results/02_cellranger_count
mkdir -p $outdir

#transcriptome to use (reference genome with only YFP and ChR2 and no PRV sequences)
transcriptome_human_noPRV=$basedir/assets/reference_genomes/refdata-lz-1266-grch38-chr2-yfp 

# run cell ranger count script for sample 7 using human reference genome with no PRV seqeunces
for d in /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/data/gi_LZ4087/LZ4087_07*_L001_R1*; do 
	#get sample name from the file name
	smp_filename=$(basename $d)
	smp=$(awk -F'_S' '{print $1}' <<< $smp_filename)

	#get sample number from the sample name
	sample_number=$(awk -F'_' '{print $2}' <<< $smp)

	#based on the sample number, assign the ouptut id, transcriptome and include_introns values
	out_id=${smp}_human_noPRV
	transcriptome_dir=$transcriptome_human_noPRV
	include_intron_val=True

	#run cellranger
	qsub $cellranger_script $out_id $fastq_dir $smp $transcriptome_dir $include_intron_val $outdir
done

################### END ###################