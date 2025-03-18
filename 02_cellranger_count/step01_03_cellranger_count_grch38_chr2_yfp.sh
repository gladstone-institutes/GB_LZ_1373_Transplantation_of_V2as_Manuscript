#!/bin/bash

# path to the cell ranger script
cellranger_script=/gladstone/bioinformatics/projects/lz-1266-lana-zholudeva-deepak-srivastava-scrnaseq-hg38-jan-2022/scripts/GB-LZ-1266/src/bash/cellranger_count_sc_sn_rna_seq.sh

#make the results output directory
outdir=/gladstone/bioinformatics/projects/lz-1266-lana-zholudeva-deepak-srivastava-scrnaseq-hg38-jan-2022/results/01_cellranger_count_grch38_chr2_yfp
mkdir -p $outdir

#transcriptome to use
transcriptome_dir=/gladstone/bioinformatics/projects/lz-1266-lana-zholudeva-deepak-srivastava-scrnaseq-hg38-jan-2022/assets/reference_genomes/refdata-lz-1266-grch38-chr2-yfp

# qsub cell ranger count script for all 5 samples 
for d in /gladstone/bioinformatics/projects/lz-1266-lana-zholudeva-deepak-srivastava-scrnaseq-hg38-jan-2022/data/gi-LZ3672/*_R1*; do
	smp_filename=$(basename $d)
	smp=$(awk -F'_S' '{print $1}' <<< $smp_filename)
	qsub $cellranger_script $smp $transcriptome_dir $outdir
done

################### END ###################
