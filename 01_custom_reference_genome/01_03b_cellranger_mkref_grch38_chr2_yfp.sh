#!/bin/bash

#genome metadata
genome="refdata-lz-1373-grch38-chr2-yfp"
version="2022-A"

#paths
build_dir=/gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/results/01_custom_reference_genome/reference_genomes_intermediate_files/refdata-lz-1373-grch38-chr2-yfp-2022-A_build
output_dir=/gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/results/01_custom_reference_genome/reference_genomes
mkdir -p "$output_dir"

#input files
fasta_custom="${build_dir}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.custom"
gtf_custom="${build_dir}/gencode.v32.primary_assembly.annotation.gtf.custom"

cd $output_dir
cellranger mkref --ref-version="$version" --genome="$genome" --fasta="$fasta_custom" --genes="$gtf_custom"


################### END ###################