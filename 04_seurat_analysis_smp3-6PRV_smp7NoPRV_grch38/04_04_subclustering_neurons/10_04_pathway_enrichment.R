#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-LZ-1373
## Authors: Ayushi Agrawal
##
## Script Goal: Perform pathway enrichment
##
## Usage example:
## Rscript 10_04_pathway_enrichment.R \
## --input clustering_cell_type/20PC_0.2res \ # Clustering results directory
## --output output_dir \                      # Output directory
## --p_val_cutoff 0.01 \                      # DE adjusted p value cutoff 
## --fold_change_cutoff 1 \                   # DE log fold change cutoff
## --pathway_db Hs_20240615.RData \           # RData file of pathway DBs
## --species 'human'                          # Species ('human' or 'mouse')
##
## Run "10_04_pathway_enrichment.R --help" for more information
###############################################################################

# Get input arguments -----------------------------------------------------
library(optparse)
option_list <- list(
  make_option(c("-i", "--input"),
              action = "store", default = NA, type = "character",
              help = "Input directory with clustering results (required)"
  ),
  make_option(c("-o", "--output"),
              action = "store", default = NA, type = "character",
              help = "Output directory, will create if it doesn't exist (required)"
  ),
  make_option("--minGSSize",
              action = "store", default = 5, type = "numeric",
              help = "Min set size for pathway enrichment (optional)"
  ),
  make_option("--maxGSSize",
              action = "store", default = 500, type = "numeric",
              help = "Max set size for pathway enrichment (optional)"
  ),
  make_option("--pathway_db",
              action = "store", default = NA, type = "character",
              help = "RData file that contains pathway DBs in the same format used by
    Interactive Enrichment Analysis (required)"
  ),
  make_option("--p_val_cutoff",
              action = "store", default = 0.05, type = "numeric",
              help = "DE adjusted p value cutoff for inclusion in ORA (optional)"
  ),
  make_option("--fold_change_cutoff",
              action = "store", default = 0.25, type = "numeric",
              help = "DE log fold change cutoff for inclusion in ORA (optional)"
  ),
  make_option("--species",
              action = "store", default = "human", type = "character",
              help = "Species - supported options are 'human' or 'mouse',
              default is 'human'"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))


# Check if required args are provided
if (is.na(opt$input) | is.na(opt$output) | is.na(opt$pathway_db) | is.na(opt$species)) {
  stop("Missing one or more required arguments")
}

# Create the results folders
if (!(dir.exists(opt$output))) {
  dir.create(opt$output, recursive = T)
}
# Change working directory
setwd(opt$output)




# Packages ----------------------------------------------------------------
# Load required packages
library(clusterProfiler)
library(dplyr)
library(Seurat)
library(tidyr)
library(stringr)
library(magrittr)
library(ggplot2)
library(gridExtra)
library(forcats)
library(EnhancedVolcano)
library(ggpubr)
library(gridExtra)
library(readr)
library(grid)
set.seed(42)




# Functions ---------------------------------------------------------------
# Volcano plots
plot_volcano <- function(de,
                         cell_type,
                         height = 10,
                         width = 15,
                         p_value_cutoff = opt$p_val_cutoff,
                         fold_change_cutoff = opt$fold_change_cutoff,
                         output) {
  number_of_de_genes <- de %>%
    filter(p_val_adj < p_value_cutoff & abs(avg_log2FC) > fold_change_cutoff) %>%
    nrow()
  subtitle_expr <- bquote("Log"[2] ~ "fold change cutoff " ~ .(as.character(fold_change_cutoff)) ~ "; Adjusted p-value cutoff " ~ .(as.character(p_value_cutoff)))
  plot <- EnhancedVolcano(de,
                          title = paste0(cell_type,"\nNumber of significant DE genes = ", number_of_de_genes),
                          subtitle = subtitle_expr,
                          lab = de$gene,
                          pCutoff = p_value_cutoff,
                          FCcutoff = fold_change_cutoff,
                          x = "avg_log2FC",
                          y = "p_val_adj",
                          #caption = caption_expr,
                          ylab = bquote(~-Log[10]~italic(p_val_adj)),
                          drawConnectors = TRUE,
                          widthConnectors = 0.3,
                          boxedLabels = TRUE,
                          labSize = 3.0,
                          colConnectors = "black",
                          legendLabSize = 14,
                          legendIconSize = 4.0)
  ggsave(plot = plot,
         filename = paste0( output, "/", cell_type, "_volcano.pdf"),
         width = width,
         height = height,
         units = "in")
}


# Convert ids
map_ids <- function(input, db, fromType) {
  toType.list <- c("ENTREZID", "SYMBOL")
  if (fromType == "ENTREZID") {
    toType.list <- c("SYMBOL")
  }
  
  # perform ID mapping
  output <- bitr(input$gene,
                 fromType = fromType,
                 toType = toType.list,
                 OrgDb = db
  )
  
  join_cols <- c("gene")
  names(join_cols) <- fromType
  output <- dplyr::left_join(output, input, by = join_cols)
  return(output)
}


# Function to extract the desired strings based on patterns
extract_category <- function(x) {
  if (grepl("^go_", x, ignore.case = TRUE)) {
    return("GO")
  } else if (grepl("^pfocr_", x, ignore.case = TRUE)) {
    return("PFOCR")
  } else if (grepl("^wp_", x, ignore.case = TRUE)) {
    return("WikiPathways")
  }
}


# ORA enrichment bar plot
plot_enrichment_barplot <- function(enrich_df, cell_type, output, db,
                                    height = 10,
                                    width = 18) {
  db <- extract_category(db)
  enrich_df <- enrich_df |>
    mutate(negative_log10_padj = -log10(p.adjust)) |>
    arrange(desc(negative_log10_padj)) |>
    head(n = 20)
  xlab_expr <- bquote("-Log"[10] ~ "(P"["adj"] ~ ")")
  
  p <- ggplot(enrich_df, aes(
    y = reorder(Description, negative_log10_padj),
    x = negative_log10_padj,
    fill = negative_log10_padj
  )) +
    geom_bar(stat = "identity", color = "black") +
    scale_y_discrete(position = "right") +
    scale_fill_viridis_c(option = "magma") +
    ylab(label = "") +
    xlab(label = xlab_expr) +
    labs(
      title = paste0("Cluster ",cell_type),
      subtitle = paste0("ORA Terms for ", db)
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none")
  
  ggsave(
    plot = p,
    filename = paste0( output, "/cluster_", cell_type, "_", db, "_enrichment_barplot.pdf"),
    width = width,
    height = height,
    units = "in"
  )
}

# ORA heatmap
plot_ora_heatmap <- function(enrichResult, original_df, output, db, cell_type,
                             p_value_cutoff = 0.05, top_count = 25) {
  db <- extract_category(db)
  plot_df <- enrichResult %>%
    top_n(n = top_count, wt = Count) %>%
    mutate(geneID = str_split(geneID, "/")) %>%
    unnest(cols = geneID) %>%
    dplyr::rename(gene = geneID) %>%
    inner_join(original_df %>%
                 select(-rank) %>%
                 filter(p.value < opt$p_val_cutoff & abs(fold.change) > opt$fold_change_cutoff), by = "gene") %>%
    filter(pvalue < p_value_cutoff)
  
  
  p1 <- plot_df %>%
    mutate(Description = str_trunc(Description, 60)) %>%
    complete(gene, Description) %>%
    arrange(fold.change) %>%
    ggplot(aes(
      x = fct_reorder(gene, abs(fold.change), .desc = TRUE),
      y = Description,
      fill = fold.change
    )) +
    geom_tile() +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_gradient2(
      low = "#1d80d1",
      mid = "white",
      high = "#f05316",
      midpoint = 0.0
    ) +
    labs(
      fill = "log2 Fold\nChange",
      title = paste0("ORA Terms with Top ", top_count, " Counts for ", db)
    ) +
    xlab("Genes")
  
  if (length(plot_df$gene %>% unique()) >= 40) {
    w <- 0.5 * length(plot_df$gene %>% unique())
  } else {
    w <- 25
  }
  
  ggsave(
    plot = p1,
    file.path(output, paste0( "cluster_", cell_type, "_", db, "_ORA_heatmap_top_enriched.pdf")), 
    width = w,
    limitsize = FALSE,
    height = 16,
    units = "in"
  )
}


# Run ORA
run_ora <- function(ranked_genes,
                    output,
                    fromType,
                    cell_type) {
  
  if(!dir.exists(output)){
    dir.create(output, recursive = TRUE)
  }
  
  # Identifier mapping
  ranked_genes_entrez <- map_ids(ranked_genes, org.db.name, fromType)
  
  gene_list <- ranked_genes_entrez %>%
    filter(p.value < opt$p_val_cutoff & abs(fold.change) > opt$fold_change_cutoff) %>%
    pull(ENTREZID)
  if (length(gene_list) == 0) {
    return()
  }
  
  universe <- map_ids(data.frame(gene=background_genes), org.db.name, fromType)
  
  for (this_db in db) {
    # Object from string
    database <- eval(parse(text = this_db))
    db_name <- extract_category(this_db)
    
    enrichResult <- enricher(
      gene_list,
      universe = universe$ENTREZID,
      TERM2GENE = database[, c("term", "gene")],
      TERM2NAME = database[, c("term", "name")],
      minGSSize = opt$minGSSize,
      maxGSSize = opt$maxGSSize,
      # pAdjustMethod="holm", #default is "BH"
      pvalueCutoff = 1 # to limit results
    )
    
    if (!is.null(enrichResult)) {
      enrichResult <- setReadable(enrichResult, org.db.name, keyType = "ENTREZID")
      if (!is.null(enrichResult)) {
        enrichResult <- as.data.frame(enrichResult@result)
        dir.create(file.path(output, this_db), showWarnings = F, recursive = T)
        
        write_csv(enrichResult, 
                  file.path(output, this_db, 
                            paste0("cluster_",cell_type,"_",db_name, "_ORA_results_table.csv")))
        if (nrow(enrichResult) > 0) {
          plot_ora_heatmap(
            enrichResult = enrichResult,
            original_df = ranked_genes,
            output = file.path(output, this_db),
            db = this_db,
            cell_type = cell_type
          )
          plot_enrichment_barplot(
            enrich_df = enrichResult,
            cell_type = cell_type,
            output = file.path(output, this_db),
            db = this_db
          )
        }
      }
    }
  }
  
  #-----------------------------------
  #KEGG analysis
  #------------------------------
  ekeggbp <- enrichKEGG(
    gene     = gene_list %>% subset(., !is.na(.)),
    universe = universe$ENTREZID %>% subset(., !is.na(.)),
    organism    = "hsa",
    minGSSize = 10,
    pvalueCutoff = 0.8,
    keyType = "ncbi-geneid"
  )

  #translating gene IDs to human readable symbols
  if (!is.null(ekeggbp)) {
    ekeggbp <- setReadable(ekeggbp, OrgDb = org.db.name, keyType="ENTREZID")
    if (nrow(ekeggbp) > 0) {
      #Visualize
      dir.create(file.path(output, "KEGG"), showWarnings = F, recursive = T)
      ## save images
      pdf(file.path(output, "KEGG", 
                    paste0( "cluster_", cell_type, "_KEGG_dotplot.pdf")),
          height = 8)
      print(dotplot(ekeggbp, showCategory = 20, orderBy="GeneRatio") )
      dev.off()
      
      x2 <- enrichplot::pairwise_termsim(ekeggbp)
      pdf(file.path(output, "KEGG", 
                    paste0("cluster_", cell_type, "_KEGG_erichment_map.pdf")),
          height = 12,
          width = 14)
      print(emapplot(x2, showCategory = 30, cex_label_category=0.8) )
      dev.off()
      
      #save the list of enriched pathways
      write.csv(ekeggbp,
                file = file.path(output, "KEGG", 
                                 paste0("cluster_", cell_type, "_KEGG_ORA_results_table.csv")))
    }
  }
  
}




# Run analysis ------------------------------------------------------------
if(opt$species == "human"){
  org.db.name <- "org.Hs.eg.db"
} else if(opt$species == "mouse"){
  org.db.name <- "org.Mm.eg.db"
} else {
  stop("Species must be either 'human' or 'mouse'")
}

# load pathway dbs
db <- load(opt$pathway_db)

# Read in the DE markers for all clusters 
data <- read.csv(list.files(opt$input, "marker_genes_per_cluster.csv", full.names = T))

# Read in the Seurat object
Seurat_obj <- readRDS(list.files(opt$input, ".rds", full.names = T))
Seurat_obj <- PrepSCTFindMarkers(Seurat_obj)

# create output directories
dir.create( paste0(opt$output, "/pathway_enrichment/ORA"), recursive = T, showWarnings = F)
dir.create( paste0(opt$output, "/pathway_enrichment/volcano_plots"), recursive = T, showWarnings = F)

# find the list of background genes
all_data <- GetAssayData(Seurat_obj, assay = "SCT")
background_genes <- rownames(all_data[rowSums(all_data[])>0,])
write.csv(as.data.frame(background_genes),
          file = paste0(opt$output, 
                        "/pathway_enrichment/ORA/nonzero_bakcground_gene_list.csv"),
          row.names = FALSE)


clusters <- unique(data$cluster)


for (c in 1:length(clusters)) {
  plot_volcano(data[data$cluster == clusters[c],],
               cell_type = paste0("Cluster ",clusters[c]),
               output = paste0(opt$output, "/pathway_enrichment/volcano_plots"))
  
  
  ranked_de_genes <- data[data$cluster == clusters[c],] %>%
    as.data.frame() %>%
    dplyr::select(gene, avg_log2FC, p_val_adj) %>%
    mutate(rank = sign(avg_log2FC) * -log10(p_val_adj)) %>%
    arrange(rank) %>%
    dplyr::rename(fold.change = avg_log2FC, p.value = p_val_adj)
  
  
  run_ora(
    ranked_genes = ranked_de_genes,
    output = paste0(opt$output, "/pathway_enrichment/ORA/cluster_", clusters[c]),
    fromType = "SYMBOL",
    cell_type = clusters[c]
  )
}



# save the session info
writeLines(
  capture.output(sessionInfo()),
  file.path(paste0(opt$output, "/pathway_enrichment"), "sessionInfo.txt")
)


# END --------------------------------------------------------
