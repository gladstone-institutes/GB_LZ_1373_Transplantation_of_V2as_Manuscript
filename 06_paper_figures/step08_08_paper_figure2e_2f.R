# load required libraries
library(Seurat)
library(ggplot2)

# set working directory
setwd(paste0("~/Dropbox (Gladstone)/GB-LZ-1373/results/",
             "03_scRNAseq_seurat_analysis_smp8_oldsmp/",
             "03_clustering_cell_type/15PC_0.1res/"))


# load the single-cell data
sc_dat <- readRDS("samples_8_and_oldsample_grch38_clustered_and_cell_typed.rds")


# inspect the data
table(sc_dat$orig.ident)
# SeuratProject 
# 22070

table(sc_dat$seurat_clusters)
# 0     1     2     3     4 
# 11081  5964  2472  1587   966 


# generate the dotplot
pdf(file = paste0("samples_8_and_oldsample_grch38_singlecell_15PC_0.1res_dotplot_fig2e_2f.pdf"))
print(DotPlot(sc_dat, 
              features = rev(c("PLP1", "SLC1A3", "PAX3", "PRDM8", "SOX2", 
                           "MAG", "MOG", "GFAP", "TUBB3")),
              group.by = "orig.ident",
              cols = c("#33BBEE", "#EE3377")) +
        RotatedAxis() + coord_flip())
dev.off()


## END