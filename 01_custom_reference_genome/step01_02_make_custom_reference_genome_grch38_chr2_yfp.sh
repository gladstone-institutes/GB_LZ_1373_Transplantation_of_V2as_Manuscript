#!/bin/bash
#$ -cwd
#$ -o /gladstone/bioinformatics/projects/lz-1266-lana-zholudeva-deepak-srivastava-scrnaseq-hg38-jan-2022/tmp/
#$ -e /gladstone/bioinformatics/projects/lz-1266-lana-zholudeva-deepak-srivastava-scrnaseq-hg38-jan-2022/tmp/
#$ -pe smp 2
#$ -l mem_free=50G
#$ -l scratch=50G
#$ -l h_rt=02:00:00
#$ -j yes

data_dir=/gladstone/bioinformatics/projects/lz-1266-lana-zholudeva-deepak-srivastava-scrnaseq-hg38-jan-2022/
script_dir=/gladstone/bioinformatics/projects/lz-1266-lana-zholudeva-deepak-srivastava-scrnaseq-hg38-jan-2022/scripts/GB-LZ-1266/src/bash
container_dir=/gladstone/bioinformatics/projects/lz-1266-lana-zholudeva-deepak-srivastava-scrnaseq-hg38-jan-2022/assets/containers
export SINGULARITY_BINDPATH="$data_dir"


#modify and filter source files 
$script_dir/modify_source_files_grch38.sh
echo "*********   modify done!  *************"

#add the custom marker genes (i.e. ChR2 and YFP) to the modified source files
if [ "$?" -eq 0 ]; then
    #append the marker genes sequences
    echo "*********   marker genes append starting!  *************"
    $script_dir/append_marker_gene_sequence_grch38_chr2_yfp.sh
    echo "*********   marker genes append done!  *************"
    if [ "$?" -eq 0 ]; then
    	#run the cellranger mkref to generate the custom reference genome
        singularity exec $container_dir/cellranger-v6.1.1.sif $script_dir/cellranger_mkref_grch38_chr2_yfp.sh
        echo "************* end  *************!"
    fi
fi


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
