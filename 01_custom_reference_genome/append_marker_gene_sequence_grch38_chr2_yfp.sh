#!/bin/bash

# Set up source and build directories
build="/gladstone/bioinformatics/projects/lz-1266-lana-zholudeva-deepak-srivastava-scrnaseq-hg38-jan-2022/assets/reference_genomes_intermediate_files/refdata-lz-1266-grch38-chr2-yfp-2022-A_build"
mkdir -p "$build"

# set paths for input files
source="/gladstone/bioinformatics/projects/lz-1266-lana-zholudeva-deepak-srivastava-scrnaseq-hg38-jan-2022/assets/reference_genomes_intermediate_files/grch38_reference_sources"
fasta_modified="${source}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.modified"
fasta_custom="${build}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.custom"
fasta_marker_genes="/gladstone/bioinformatics/projects/lz-1266-lana-zholudeva-deepak-srivastava-scrnaseq-hg38-jan-2022/assets/marker_gene_sequence/marker_genes.fa"
gtf_filtered="${source}/gencode.v32.primary_assembly.annotation.gtf.filtered"
gtf_custom="${build}/gencode.v32.primary_assembly.annotation.gtf.custom"
gtf_marker_genes="/gladstone/bioinformatics/projects/lz-1266-lana-zholudeva-deepak-srivastava-scrnaseq-hg38-jan-2022/assets/marker_gene_sequence/marker_genes.gtf"

#append human MAPT fasta sequence at the end of the reference fasta file
cp $fasta_modified $fasta_custom 
cat $fasta_marker_genes >> $fasta_custom
grep ">" $fasta_custom

#append human MAPT GTF file at the end of the reference GTF file
cp $gtf_filtered $gtf_custom
cat $gtf_marker_genes >> $gtf_custom
tail $gtf_custom

################### END ###################
