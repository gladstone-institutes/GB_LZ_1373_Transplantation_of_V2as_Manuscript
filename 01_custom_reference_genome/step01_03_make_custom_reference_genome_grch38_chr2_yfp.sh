#!/bin/bash
#$ -cwd
#$ -o /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/tmp/
#$ -e /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/tmp/
#$ -pe smp 2
#$ -l mem_free=50G
#$ -l scratch=50G
#$ -l h_rt=02:00:00
#$ -j yes

data_dir=/gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023
script_dir=$data_dir/scripts/GB_LZ_1373_Transplantation_of_V2as_Manuscript/01_custom_reference_genome
container_dir=$data_dir/assets
export SINGULARITY_BINDPATH="$data_dir"


#add the custom marker genes (i.e. ChR2 and YFP) to the modified source files
if [ "$?" -eq 0 ]; then
    #append the marker genes sequences
    echo "*********   marker genes append starting!  *************"
    $script_dir/01_03a_append_marker_gene_sequence_grch38_chr2_yfp.sh
    echo "*********   marker genes append done!  *************"
    if [ "$?" -eq 0 ]; then
        #run the cellranger mkref to generate the custom reference genome
        singularity exec $container_dir/cellranger-v6.1.1.sif $script_dir/01_03b_cellranger_mkref_grch38_chr2_yfp.sh
        echo "************* end  *************!"
    fi
fi


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"