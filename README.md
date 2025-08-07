# GB_LZ_1373_Transplantation_of_V2as_Manuscript
Code for analysis of single-cell/nuclei RNA-seq data in the Transplantation of V2as Manuscript (with Lana Zholudeva et al.) 

## Researchers
- Lana Zholudeva
- Ayushi Agrawal

## Experimental details
- 10X single cell RNA-seq experiment and single nucleus RNA-seq experiment
- There are 8 samples in this experiment:
    - 4 Samples: Injured female Sprague Dawley Rat Cervical 3/4 spinal cord injury with human V2a Spinal interneurons, traced with PRV614(RFP) 2 months post-transplantation
    - 1 Sample: Injured female Sprague Dawley Rat Cervical 3/4 spinal cord injury with human V2a Spinal interneurons 2 months post-transplantation. NO PRV and different batch of cells
    - 2 Samples: single cell sample of cells the same batch as the original single cell dataset. These are the cells prior to transplantation
- All samples 3-8 would be aligned to human genome with PRV and RFP and ChR2 sequences added in. 
- Metadata information:

 
## Expected folder structure
GB_LZ_1373/  
│── data  
│── results  
│── assets  
│── tmp  
│── scripts/  
│&emsp;&emsp;│── GB_LZ_1373_Transplantation_of_V2as_Manuscript  

**Note:** The analysis was performed on UCSF Wynton HPC cluster. All Dockerfiles provided in the `Docker/` directory were used to create Docker containers, which were then converted to Singularity containers (`.sif` files) for HPC compatibility. 

## Important Notes for Re-running Analysis

### Path Considerations
This analysis was performed on **Gladstone's Hive storage system** and results were shared via **Dropbox**. As a result, the scripts contain hardcoded paths that reflect this specific computing environment:

- **Shell scripts** contain paths like `/gladstone/bioinformatics/projects/lz-1373-lana-zholudeva-deepak-srivastava-scrnaseq-snrnaseq-hg38-may-2023/`
- **R scripts** contain paths like `~/Dropbox (Gladstone)/GB-LZ-1373/results/`
- **Container files** (`.sif` files) are referenced but not included in this repository

### To Re-run the Analysis
If you want to re-run this analysis on your own system, you will need to:

1. **Update all hardcoded paths** in the scripts to match your local directory structure
2. **Build the required Singularity containers** using the Dockerfiles provided in the `Docker/` directory:
   - `r_seurat_4.1_harmony_gb_lz_1266.sif` (from `Dockerfile-r_seurat_4.1_harmony_gb_lz_1266`)
   - `seurat-4-3_soupx-1-6-1_doubletfinder_67fb8b5.sif` (from `Dockerfile-seurat-4-3_soupx-1-6-1_doubletfinder`)
   - `cellranger-v6.1.1.sif` (from `Dockerfile-cellranger`)
3. **Provide your own input data** (FASTQ files, reference genomes, etc.)

### Missing Components
The following components are referenced in the scripts but not included in this repository:
- **Singularity container files** (`.sif` files) - however, the Dockerfiles to build them are provided in the `Docker/` directory
- **Input data files** (FASTQ files, processed `.rds` files)

This repository serves as a **reference implementation** of the analysis pipeline used in the manuscript, but requires adaptation for use on different computing environments.

## References
1. https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
2. https://satijalab.org/seurat/articles/pbmc3k_tutorial.html


</br>
