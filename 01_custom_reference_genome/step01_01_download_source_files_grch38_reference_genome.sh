#this script should be run on the data transfer (dt) node on wynton
#!/bin/bash

# Create the reference_sources folder for the source files
source="/gladstone/bioinformatics/projects/lz-1266-lana-zholudeva-deepak-srivastava-scrnaseq-hg38-jan-2022/assets/reference_genomes_intermediate_files/grch38_reference_sources"
mkdir -p "$source"

#specify the source files to be downloaded
fasta_url="http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
fasta_download="${source}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz"
gtf_download="${source}/gencode.v32.primary_assembly.annotation.gtf.gz"

# Download source files if they do not exist in reference_sources folder
cd $source
if [ ! -f "$fasta_download" ]; then
    curl -sSO "$fasta_url"
    if [ "$?" -eq 0 ]; then
        echo "Fasta file download completed!"
    else
        echo "ERROR: Fasta file download failed."
    fi
else
	echo "Fasta file already exists."
fi
if [ ! -f "$gtf_download" ]; then
    curl -sSO "$gtf_url"
    if [ "$?" -eq 0 ]; then
        echo "GTF file download completed!"
    else
        echo "ERROR: GTF file download failed."
    fi
else
	echo "GTF file already exists."
fi


################### END ###################
