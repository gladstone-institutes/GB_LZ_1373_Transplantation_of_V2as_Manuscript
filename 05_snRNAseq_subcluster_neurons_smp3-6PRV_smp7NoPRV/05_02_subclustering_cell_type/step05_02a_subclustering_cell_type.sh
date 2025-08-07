#!/usr/bin/env bash
#$ -cwd
#$ -o /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/tmp/
#$ -e /gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/tmp/
#$ -pe smp 2
#$ -l mem_free=60G
#$ -l scratch=100G
#$ -l h_rt=50:00:00
#$ -j yes


data_dir=/gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023
script_dir=$data_dir/scripts/GB_LZ_1373_Transplantation_of_V2as_Manuscript/05_snRNAseq_subcluster_neurons_smp3-6PRV_smp7NoPRV
results_dir=$data_dir/results/05_snRNAseq_subcluster_neurons_smp3-6PRV_smp7NoPRV
container_dir=$data_dir/assets
export APPTAINER_BINDPATH="$data_dir"


# Define the resolutions to loop over
resolutions=(0.02 0.04 0.06 0.08 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 1.10 1.20)


# 15 PCs
# Loop through each resolution and run the neuron subclustering
for res in "${resolutions[@]}"; do
    singularity exec $container_dir/seurat-4-3_soupx-1-6-1_doubletfinder_67fb8b5.sif \
    Rscript $script_dir/05_02_subclustering_cell_type/05_02a_subclustering_cell_type.R \
    --input $results_dir/01_subset_and_visualize/neurons_smp3-7_grch38_noCB_noDF_NoPRV_data_sct_pca_umap.rds \
    --output $results_dir/02_subclustering_cell_type/15PC_${res}res_rmBatchEffect \
    --output_prefix "neurons_smp3-7_grch38_noCB_noDF_NoPRV_15PC_${res}res_rmBatchEffect" \
    --tissue "Brain" \
    --ndim 15 \
    --resolution "$res" \
    --batch_var "library_prep_batch"
done


# 30 PCs
# Loop through each resolution and run the neuron subclustering
for res in "${resolutions[@]}"; do
    singularity exec $container_dir/seurat-4-3_soupx-1-6-1_doubletfinder_67fb8b5.sif \
    Rscript $script_dir/05_02_subclustering_cell_type/05_02a_subclustering_cell_type.R \
    --input $results_dir/01_subset_and_visualize/neurons_smp3-7_grch38_noCB_noDF_NoPRV_data_sct_pca_umap.rds \
    --output $results_dir/02_subclustering_cell_type/30PC_${res}res_rmBatchEffect \
    --output_prefix "neurons_smp3-7_grch38_noCB_noDF_NoPRV_30PC_${res}res_rmBatchEffect" \
    --tissue "Brain" \
    --ndim 30 \
    --resolution "$res" \
    --batch_var "library_prep_batch"
done


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

######################## END ########################