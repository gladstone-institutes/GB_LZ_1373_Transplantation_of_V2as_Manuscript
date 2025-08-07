#!/usr/bin/env Rscript

###############################################################################
## Project ID: GB-LZ-1373
## Authors: Reuben Thomas, Natalie Elphick, Ayushi Agrawal
##
## Script Goal: Perform pseudo-bulked differential gene expression (one vs other clusters)
##
## Usage example:
## Rscript step05_03_pseudobulkDE.R
###############################################################################

# ------------------------------------------------------------------------------
# Script Options (previously defined via optparse)
#
# opt$input              : Input Seurat object in RDS format (required)
# opt$output             : Output directory where results will be saved. Will be created if it doesn't exist (required)
# opt$output_prefix      : Prefix for output files (optional, helps organize outputs)
# opt$cluster_of_interest: Cluster ID to compare against all other clusters in differential expression (required)
# opt$predictor          : Comma-separated predictor variables for the design formula (required)
# opt$cell_annotation    : Cell-type or cluster annotation column in metadata used for pseudo-bulk aggregation 
#                          (default = "seurat_clusters")
# opt$sampleid           : Metadata column name that identifies biological replicates
# opt$minCells           : Minimum median number of cells per replicate in a cluster to keep it (default = 10)
# opt$minGSSize          : Minimum gene set size to consider in pathway enrichment analysis (default = 5)
# opt$maxGSSize          : Maximum gene set size to consider in pathway enrichment analysis (default = 500)
# opt$pathway_db         : Path to RData file containing pathway databases in the format used by Interactive Enrichment Analysis (required)
# opt$p_val_cutoff       : Adjusted p-value threshold for selecting differentially expressed genes for ORA (default = 0.05)
# opt$fold_change_cutoff : Log fold change threshold for selecting genes in ORA (default = 0.25)
# opt$species            : Species for pathway enrichment ('human' or 'mouse'; default = 'human')
# ------------------------------------------------------------------------------

# Define input arguments -------------------------------------------------------
base_dir <- "~/Dropbox (Gladstone)/GB-LZ-1373/results/05_snRNAseq_subcluster_neurons_smp3-6PRV_smp7NoPRV"
opt <- data.frame(input = file.path(base_dir,
                                    "02_subclustering_cell_type/30PC_0.30res_rmBatchEffect",
                                    paste0("split_subcluster6/neurons_smp3-7_grch38_noCB_noDF_NoPRV_30PC_0.30res_",
                                           "rmBatchEffect_split_subcluster6_clustered_and_cell_typed_split.rds")), 
                  output = file.path(base_dir,
                                     "03_pseudobulkDE_neuron_subclustering_NoPRV"), 
                  output_prefix = "neurons", 
                  cell_annotation = "modified_clusters",
                  species = "human",
                  sampleid = "sample_name",
                  predictor = "orig.ident",
                  minCells = 10,
                  minGSSize = 5,
                  maxGSSize = 500,
                  pathway_db = "~/Dropbox (Gladstone)/github/GB_LZ_1373_Transplantation_of_V2as_Manuscript/reference/Hs_20250415.RData",
                  p_val_cutoff = 0.01,
                  fold_change_cutoff = 2
)

# Check if required args are provided
if (is.na(opt$input) | is.na(opt$output) | is.na(opt$predictor) | 
    is.na(opt$sampleid) | is.na(opt$pathway_db) | is.na(opt$species)) {
  stop("Missing one or more required arguments")
}

# create output directory
if (!(dir.exists(opt$output))) {
  dir.create(opt$output, recursive = T, showWarnings = F)
}


# Setup ------------------------------------------------------------------------
# Load required packages
library(clusterProfiler)
library(Seurat)
library(HGNChelper)
library(muscat)
library(SummarizedExperiment)
library(emmeans)
library(gridExtra)
library(EnhancedVolcano)
library(edgeR)
library(ggpubr)
library(grid)
library(openxlsx)
library(tidyverse)
library(magrittr)

# set seed
set.seed(1234)

# Set default colors to Wong B (https://doi.org/10.1038/nmeth.1618)
# colorblindness friendly colors
custom_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
  "#666666", "#AD7700", "#1C91D4", "#007756", "#D5C711", "#005685", "#A04700",
  "#B14380", "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71", "#06A5FF",
  "#FF8320", "#D99BBD", "#8C8C8C", "#FFCB57", "#9AD2F2", "#2CFFC6", "#F6EF8E",
  "#38B7FF", "#FF9B4D", "#E0AFCA", "#A3A3A3", "#8A5F00", "#1674A9", "#005F45",
  "#AA9F0D", "#00446B", "#803800", "#8D3666", "#3D3D3D"
)
options(
  ggplot2.continuous.colour = list(custom_colors),
  ggplot2.continuous.fill = list(custom_colors),
  ggplot2.discrete.colour = list(custom_colors),
  ggplot2.discrete.fill = list(custom_colors)
)


# Functions --------------------------------------------------------------------
# refine the metadata levels in a Seurat object
refine_metadata_levels <- function(seurat_data) {
  for (i in base::colnames(seurat_data@meta.data)) {
    if (base::is.factor(seurat_data@meta.data[[i]])) {
      base::print(base::paste("Re-evaluating levels for a factor column", i))
      base::print(
        base::paste(
          "before:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse = ", ")
        )
      )
      seurat_data@meta.data[[i]] <- base::droplevels(seurat_data@meta.data[[i]]) # need to drop levels of the removed values
      base::print(
        base::paste(
          "after:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse = ", ")
        )
      )
    }
  }
  return(seurat_data)
}

# Function to extract the desired strings based on patterns
extract_category <- function(x) {
  if (grepl("^go_", x, ignore.case = TRUE)) {
    return("GO")
  } else if (grepl("^pfocr_", x, ignore.case = TRUE)) {
    return("PFOCR")
  } else if (grepl("^wp_", x, ignore.case = TRUE)) {
    return("WikiPathways")
  } else if (grepl("^KEGG", x, ignore.case = TRUE)) {
    return("KEGG")
  }
}

multipage_plot <- function(plot_list, per_page, filename, ncol = 1) {
  # Create a PDF file to save the plots
  pdf(file = filename, height = 10, width = 14)
  
  # Split the plots into groups of per_page per page
  num_plots <- length(plot_list)
  num_pages <- ceiling(num_plots / per_page)
  plot_indices <- split(
    seq_len(num_plots),
    rep(seq_len(num_pages), each = per_page)
  ) |>
    suppressWarnings()
  
  # Plot the plots on each page
  for (i in seq_len(num_pages)) {
    if (length(plot_indices[[i]]) < per_page) {
      # Fill in the rest of the last page with blank plots
      this_index <- max(plot_indices[[i]])
      for (j in seq_len(per_page - length(plot_indices[[i]]))) {
        this_index <- this_index + 1
        plot_indices[[this_index]] <- ggplot() +
          geom_blank()
        plot_indices[[i]] <- c(plot_indices[[i]], this_index)
      }
    }
    plots <- ggarrange(
      plotlist = plot_list[plot_indices[[i]]],
      ncol = ncol, nrow = ceiling(length(plot_indices[[i]]) / ncol)
    )
    plots_grob <- ggplotGrob(plots)
    grid.newpage()
    grid.draw(plots_grob)
  }
  # Close the PDF file
  dev.off()
}

# Volcano plots
plot_volcano <- function(de,
                         cell_type,
                         comparison,
                         height = 10,
                         width = 15,
                         p_value_cutoff = 0.05,
                         fold_change_cutoff = 0.25,
                         output) {
  number_of_de_genes <- de %>%
    filter(FDR < p_value_cutoff & abs(logFC) > fold_change_cutoff) %>%
    nrow()
  
  filename <- paste0(cell_type, "_", comparison, ".pdf")
  subtitle_expr <- paste0("\nFDR cutoff = ", p_value_cutoff, 
                          "; log2 fold change cutoff = ", fold_change_cutoff,
                          "\nNumber of DE genes = ", number_of_de_genes)
  
  p_value_cutoff_new <- min(de$PValue[de$FDR >= p_value_cutoff])
  
  plot <- EnhancedVolcano(de,
                          title = paste0(cell_type, "\n", comparison),
                          subtitle = subtitle_expr,
                          lab = de$gene,
                          pCutoff = p_value_cutoff_new,
                          FCcutoff = fold_change_cutoff,
                          x = 'logFC',
                          y = 'PValue',
                          ylab = bquote(~-Log[10]~italic(PValue)),
                          drawConnectors = TRUE,
                          widthConnectors = 0.3,
                          boxedLabels = TRUE,
                          labSize = 3.0,
                          colConnectors = "black",
                          legendLabSize = 10,
                          legendIconSize = 3.0,
                          legendLabels = c("NS", 
                                           bquote("Log"[2]~"FC > " ~ .(fold_change_cutoff)),
                                           bquote("FDR < " ~ .(p_value_cutoff)),
                                           bquote("FDR < " ~ .(p_value_cutoff) ~ " & Log"[2]~"FC > " ~ .(fold_change_cutoff)))
  )
  ggsave(
    plot = plot,
    filename = file.path(output, paste0("volcano_plot_", filename)),
    width = width,
    height = height,
    units = "in"
  )
}

# Box plot for genes with lowest p-values
plot_boxplots_top_genes <- function(dge, 
                                    topRes, 
                                    top_n = 10, 
                                    group_var = "cluster_compare", 
                                    output = NULL, 
                                    filename = NULL) {
  # Get top N genes by p-value
  top_genes <- topRes %>%
    arrange(PValue) %>%
    head(top_n) %>%
    pull(gene)
  
  # Extract logCPM values
  logCPM <- cpm(dge, log = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    filter(gene %in% top_genes) %>%
    column_to_rownames("gene") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample")
  
  # Extract sample metadata
  group_info <- dge$samples %>%
    rownames_to_column("sample") %>%
    dplyr::select(sample, sample_name, group = all_of(group_var))
  
  # Merge logCPM and group info
  plot_df <- left_join(logCPM, group_info, by = "sample")
  
  # Reshape to long format
  plot_df_long <- plot_df %>%
    pivot_longer(cols = all_of(top_genes), names_to = "gene", values_to = "logCPM")
  
  # Plot
  p <- ggplot(plot_df_long, aes(x = group, y = logCPM)) +
    geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.6, linewidth = 0.3) +
    geom_jitter(aes(color = sample_name), width = 0.2, size = 1.5, alpha = 0.8) +
    facet_wrap(~gene, scales = "free_y", ncol = 4) +
    labs(x = group_var, y = expression(Log[2]~CPM), title = paste("Top", top_n, "Genes with Lowest p-values")) +
    theme_bw()
  
  # Save or show plot
  if (!is.null(output) && !is.null(filename)) {
    dir.create(output, showWarnings = FALSE, recursive = TRUE)
    ggsave(filename = file.path(output, filename), plot = p, width = 10, height = 10)
  } else {
    print(p)
  }
  
  
  # Create list of ggplots
  plot_list <- lapply(unique(plot_df_long$gene), function(gene_i) {
    ggplot(filter(plot_df_long, gene == gene_i), aes(x = group, y = logCPM)) +
      geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.6, linewidth = 0.3) +
      geom_jitter(aes(color = sample_name), width = 0.2, size = 1.5, alpha = 0.8) +
      labs(
        title = gene_i,
        x = group_var,
        y = expression(Log[2]~CPM)
      ) +
      theme_bw()
  })
  
  # Create output directory if needed
  if (!is.null(output)) dir.create(output, showWarnings = FALSE, recursive = TRUE)
  
  # Use your multipage_plot function to save
  multipage_plot(plot_list, per_page = 6,
                 filename = file.path(output, filename), ncol = 2)
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

# Run ORA
run_ora <- function(gene_list,
                    output,
                    output_prefix,
                    fromType
) {
  
  # Identifier mapping
  ranked_genes_entrez <- map_ids(gene_list, org.db.name, fromType)
  
  gene_list_entrez <- ranked_genes_entrez %>%
    filter(p.value < opt$p_val_cutoff & abs(fold.change) > opt$fold_change_cutoff) %>%
    pull(ENTREZID)
  
  if (length(gene_list) == 0) {
    return()
  }
  
  # create output directory
  if(!dir.exists(output)){
    dir.create(output, recursive = TRUE)
  }
  
  universe <- map_ids(data.frame(gene = background_genes), org.db.name, fromType)
  
  # Backup existing ggplot2 options
  old_options <- list(
    ggplot2.continuous.colour = getOption("ggplot2.continuous.colour"),
    ggplot2.continuous.fill = getOption("ggplot2.continuous.fill"),
    ggplot2.discrete.colour = getOption("ggplot2.discrete.colour"),
    ggplot2.discrete.fill = getOption("ggplot2.discrete.fill")
  )
  
  # Temporarily unset ggplot2 options for dotplot
  options(
    ggplot2.continuous.colour = NULL,
    ggplot2.continuous.fill = NULL,
    ggplot2.discrete.colour = NULL,
    ggplot2.discrete.fill = NULL
  )
  
  for (this_db in db) {
    # Object from string
    database <- eval(parse(text = this_db))
    db_name <- paste0("ORA_",extract_category(this_db))
    
    enrichResult <- enricher(
      gene_list_entrez,
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
      
      if (nrow(enrichResult) > 0) {
        dir.create(file.path(output, this_db), showWarnings = F, recursive = T)
        
        write_csv(as.data.frame(enrichResult@result), 
                  file.path(output, this_db, 
                            paste0(db_name, 
                                   "_",
                                   output_prefix, 
                                   ".csv")))
        
        #trim descriptions to 80 characters
        enrichResult@result <- dplyr::mutate(enrichResult@result, 
                                             Description = str_trunc(Description, 80))
        
        p <- barplot(enrichResult, 
                     showCategory = 20,
                     label_format=50) +
          labs(
            title = output_prefix,
            subtitle = paste0("\nORA Terms for ", this_db)
          )
        ggplot2::ggsave(p, file = file.path(output, 
                                            this_db,
                                            paste0(db_name, 
                                                   "_enrichment_barplot_", 
                                                   output_prefix, 
                                                   ".pdf")), 
                        width = 2400, height = 2600, units = "px", device='pdf')
        
        p <- enrichplot::dotplot(enrichResult, 
                                 showCategory = 20,
                                 label_format=50) +
          labs(
            title = output_prefix,
            subtitle = paste0("\nORA Terms for ", this_db)
          )
        ggplot2::ggsave(p, file = file.path(output,  
                                            this_db,
                                            paste0(db_name, 
                                                   "_enrichment_dotplot_", 
                                                   output_prefix, 
                                                   ".pdf")), 
                        width = 2400, height = 2600, units = "px", device='pdf')
        
      }
    }
  }
  
  #-----------------------------------
  #KEGG analysis
  #-----------------------------------
  print("Starting KEGG...")
  ekeggbp <- enrichKEGG(
    gene     = gene_list_entrez %>% subset(., !is.na(.)),
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
      
      #save the list of enriched pathways
      write.csv(ekeggbp,
                file = file.path(output, "KEGG", 
                                 paste0( "ORA_KEGG_results_table_", 
                                         output_prefix, 
                                         ".csv")))
      
      ## save images
      #trim descriptions to 80 characters
      ekeggbp@result <- dplyr::mutate(ekeggbp@result, 
                                      Description = str_trunc(Description, 80))
      
      p <- barplot(ekeggbp, 
                   showCategory = 20,
                   label_format=50) +
        labs(
          title = output_prefix,
          subtitle = "\nORA Terms for KEGG"
        )
      ggplot2::ggsave(p, file = file.path(output, 
                                          "KEGG",
                                          paste0("ORA_KEGG_enrichment_barplot_", 
                                                 output_prefix, 
                                                 ".pdf")), 
                      width = 2400, height = 2600, units = "px", device='pdf')
      
      p <- enrichplot::dotplot(ekeggbp, 
                               showCategory = 20,
                               label_format=50) +
        labs(
          title = output_prefix,
          subtitle = "\nORA Terms for KEGG"
        )
      ggplot2::ggsave(p, 
                      file = file.path(output, 
                                       "KEGG",
                                       paste0("ORA_KEGG_enrichment_dotplot_", 
                                              output_prefix, 
                                              ".pdf")), 
                      width = 2400, height = 2600, units = "px", device='pdf')
      
    }
  }
  
  # Restore the original ggplot2 options
  options(old_options)
  
}

# Run GSEA
run_gsea <- function(ranked_genes, 
                     output, 
                     fromType,
                     output_prefix, 
                     seed = 1) {
  set.seed(seed)
  minGSSize <- opt$minGSSize
  maxGSSize <- opt$maxGSSize
  
  # Identifier mapping
  ranked_genes_entrez <- map_ids(ranked_genes, org.db.name, fromType)
  
  # Record unmapped rows
  ranked_genes_unmapped <- ranked_genes %>%
    dplyr::filter(!gene %in% ranked_genes_entrez[[fromType]]) %>%
    dplyr::arrange(desc(rank))
  
  # Resolve duplicates (keep ENTREZID with largest abs(rank))
  ranked_genes_entrez_dedup <- ranked_genes_entrez %>%
    dplyr::mutate(absrank = abs(rank)) %>%
    dplyr::arrange(desc(absrank)) %>%
    dplyr::distinct(ENTREZID, .keep_all = T) %>%
    dplyr::select(-absrank) %>%
    dplyr::arrange(desc(rank))
  
  ranked_genes_entrez_l <- ranked_genes_entrez_dedup$rank
  names(ranked_genes_entrez_l) <- ranked_genes_entrez_dedup$ENTREZID
  ranked_genes_entrez_l <- na.omit(ranked_genes_entrez_l)
  
  geneList <- sort(ranked_genes_entrez_l, decreasing = T)
  
  # Backup existing ggplot2 options
  old_options <- list(
    ggplot2.continuous.colour = getOption("ggplot2.continuous.colour"),
    ggplot2.continuous.fill = getOption("ggplot2.continuous.fill"),
    ggplot2.discrete.colour = getOption("ggplot2.discrete.colour"),
    ggplot2.discrete.fill = getOption("ggplot2.discrete.fill")
  )
  
  # Temporarily unset ggplot2 options for dotplot
  options(
    ggplot2.continuous.colour = NULL,
    ggplot2.continuous.fill = NULL,
    ggplot2.discrete.colour = NULL,
    ggplot2.discrete.fill = NULL
  )
  
  for (this_db in db) {
    print(this_db)
    # Object from string
    database <- eval(parse(text = this_db))
    
    # Perform GSEA
    gseaResult <- clusterProfiler::GSEA(
      geneList,
      TERM2GENE = database[, c("term", "gene")],
      TERM2NAME = database[, c("term", "name")],
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      # pAdjustMethod="holm", #default is "BH"
      pvalueCutoff = 0.99, # to limit results
      verbose = F
    )
    
    
    if (!is.null(gseaResult) && nrow(gseaResult) > 0) {
      gseaResult <- setReadable(gseaResult,
                                OrgDb = org.db.name,
                                keyType = "ENTREZID"
      )
      
      dir.create(file.path(output, this_db), showWarnings = F, recursive = T)
      
      write_csv(as.data.frame(gseaResult@result), 
                file.path(output, 
                          this_db, 
                          paste0("GSEA_",
                                 extract_category(this_db), 
                                 "_results_table_",
                                 output_prefix, 
                                 ".csv")))
      
      #trim descriptions to 80 characters
      gseaResult@result <- dplyr::mutate(gseaResult@result, 
                                         Description = str_trunc(Description, 80))
      
      p <- enrichplot::dotplot(gseaResult, 
                               showCategory = 20,
                               label_format=50) +
        labs(
          title = output_prefix,
          subtitle = paste0("\nORA Terms for ", this_db)
        )
      ggplot2::ggsave(p, file = file.path(output,  
                                          this_db,
                                          paste0("GSEA_",
                                                 extract_category(this_db), 
                                                 "_enrichment_dotplot_", 
                                                 output_prefix, 
                                                 ".pdf")), 
                      width = 2400, height = 2600, units = "px", device='pdf')
    }
  }
  
  #-----------------------------------
  #KEGG analysis
  #------------------------------
  ekeggbp <- gseKEGG(geneList = geneList,
                     organism = 'hsa',
                     minGSSize = minGSSize,
                     maxGSSize = maxGSSize,
                     pvalueCutoff = 0.99,
                     verbose = FALSE)
  if (!is.null(ekeggbp)) {
    gseaResult <- setReadable(ekeggbp,
                              OrgDb = org.db.name,
                              keyType = "ENTREZID"
    )
  }
  
  if (nrow(gseaResult) > 0) {
    #Visualize
    dir.create(file.path(output, "KEGG"), showWarnings = F, recursive = T)
    
    write_csv(gseaResult@result, 
              file.path(output, 
                        "KEGG", 
                        paste0("GSEA_KEGG_results_table_",
                               output_prefix, 
                               ".csv")))
    ## save images
    #trim descriptions to 80 characters
    gseaResult@result <- dplyr::mutate(gseaResult@result, 
                                       Description = str_trunc(Description, 80))
    
    p <- enrichplot::dotplot(gseaResult, 
                             showCategory = 20,
                             label_format=50) +
      labs(
        title = output_prefix,
        subtitle = "\nORA Terms for KEGG"
      )
    ggplot2::ggsave(p, 
                    file = file.path(output, 
                                     "KEGG",
                                     paste0("GSEA_KEGG_enrichment_dotplot_", 
                                            output_prefix, 
                                            ".pdf")), 
                    width = 2400, height = 2600, units = "px", device='pdf')
    
  }
  
  # Restore the original ggplot2 options
  options(old_options)
  
}


# Process inputs ---------------------------------------------------------------
if(opt$species == "human"){
  org.db.name <- "org.Hs.eg.db"
} else if(opt$species == "mouse"){
  org.db.name <- "org.Mm.eg.db"
} else {
  stop("Species must be either 'human' or 'mouse'")
}

# load pathway dbs
db <- load(opt$pathway_db)

# Read in the Seurat object
data <- readRDS(opt$input)

# sample counts per cluster
smp_cluster_counts <- unique(data@meta.data %>%
                               group_by(sample_name) %>%
                               mutate(total_neurons_per_sample = 
                                        n()) %>%
                               group_by(seurat_clusters, .add = TRUE) %>%
                               mutate(neurons_per_sample_in_cluster = 
                                        n()) %>%
                               select(sample_name,
                                      seurat_clusters,
                                      total_neurons_per_sample,
                                      neurons_per_sample_in_cluster))
write.csv(smp_cluster_counts, 
          file = file.path(opt$output,
                           paste0(opt$output_prefix,
                                  "_counts_per_sample_per_subcluster.csv")),
          row.names = FALSE)

# process input data
data <- subset(data, modified_clusters %in% c("6_1","6_2"))
data <- subset(data, sample_name != "LZ4087_07" )
data <- refine_metadata_levels(data) # Fix the levels

# Predictors
predictors <- opt$predictor
predictors <- strsplit(predictors, ",")[[1]] %>% trimws()

# Sample identifier in metadata
sampleid <- opt$sampleid %>% trimws()

# Cluster annotation column
clusters <- opt$cell_annotation
data@meta.data[, "clusters_compare"] <- data@meta.data[, clusters]

# Store the meta-data for each cell in the PhenoData object
PhenoData <- data@meta.data
PhenoData <- PhenoData[,!grepl("^DF\\.classifications_", colnames(PhenoData))]
PhenoData <- PhenoData[,!grepl("^pANN_", colnames(PhenoData))]

# Check if sampleid is in colnames of PhenoData
if (!(sampleid %in% colnames(PhenoData))) {
  stop("Columns of the meta.data slot of the Seurat object should include sampleid")
}

# get the cluster levels
cluster_levels <- PhenoData[[clusters]] %>%
  unique() %>%
  sort()


# create output directories ----------------------------------------------------
# Define the output directories in a named list
dirs <- list(
  de_outputs = file.path(opt$output, "gene_expression_associations"),
  ora_outputs = file.path(opt$output, "ORA"),
  gsea_outputs = file.path(opt$output, "GSEA")
)

# Create the directories recursively
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# Assign the paths to variables in the global environment
invisible(list2env(dirs, envir = .GlobalEnv))



# Pseudo bulk for DE analysis --------------------------------------------------
# find the list of background genes for ORA analysis later
data <- PrepSCTFindMarkers(data)
all_data <- GetAssayData(data, assay = "SCT")
background_genes <- rownames(all_data[rowSums(all_data[])>0,])
write.csv(as.data.frame(background_genes),
          file = file.path(opt$output, 
                           "ORA",
                           paste0(opt$output_prefix, "_nonzero_background_gene_list.csv")),
          row.names = FALSE)

# Create SingleCellExperiment object
sce <- SummarizedExperiment(assays = list(
  counts = data[["RNA"]]$counts
), colData = PhenoData)

sce <- as(sce, "SingleCellExperiment")
rm(data)


# Prep this object for subsequent aggregation analyses
sce <- prepSCE(sce,
               kid = clusters, # sub-population assignments
               gid = predictors[1], # group IDs
               sid = sampleid, # sample IDs
               drop = FALSE
)

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids
names(sids) <- sids

# keep those clusters with a median of at least 10 cells across all individual replicates
toKeep <- table(sce$cluster_id, sce$sample_id) %>%
  t() %>%
  apply(., 2, function(x) median(x)) %>%
  subset(., is_greater_than(., opt$minCells)) %>%
  names()

if (!("6_2" %in% toKeep)) {
  stop(print(paste0(
    "Pseudo-bulk differential expression analyses on subcluster ",
    "6_2",
    " not performed. The median number of cells across all biological replicates is less than 10 for this subcluster."
  )))
}

# Aggregate counts across cells for each sample_id within each cluster (cluster_id)
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("clusters_compare", "sample_id"))


# Perform DE analysis 
# Prepare for edgeR DE analysis
counts1 <- pb@assays@data$`6_2`
colnames(counts1) <- paste0(colnames(counts1), "_cluster_6_2")
counts2 <- pb@assays@data$`6_1`
colnames(counts2) <- paste0(colnames(counts2), "_cluster_6_1")

sampleinfo1 <- cbind(colData(pb), 
                     cluster_compare = "6_2", 
                     sample_name = rownames(colData(pb)))
rownames(sampleinfo1) <- paste0(colnames(counts1))
sampleinfo2 <- cbind(colData(pb), 
                     cluster_compare = "6_1", 
                     sample_name = rownames(colData(pb)))
rownames(sampleinfo2) <- paste0(colnames(counts2))

# Strip suffixes to get original sample names
sample_names1 <- gsub("_cluster_6_2$", "", colnames(counts1))
sample_names2 <- gsub("_cluster_6_1$", "", colnames(counts2))

# Identify samples with zero counts in either counts1 or counts2
zero_samples <- unique(c(
  sample_names1[colSums(counts1) == 0],
  sample_names2[colSums(counts2) == 0]
))
if(length(zero_samples) > 0){
  print(paste0("** Remove samples with zero total counts: ", zero_samples))
}

# Remove zero-count samples from both counts1 and counts2
counts1 <- counts1[, !sample_names1 %in% zero_samples]
counts2 <- counts2[, !sample_names2 %in% zero_samples]

# Update sample information to match filtered counts
sampleinfo1 <- sampleinfo1[!rownames(sampleinfo1) %in% paste0(zero_samples, "_cluster_6_2"), ]
sampleinfo2 <- sampleinfo2[!rownames(sampleinfo2) %in% paste0(zero_samples, "_cluster_6_1"), ]


# Combine counts and sample information for DGEList creation
dge <- DGEList(cbind(counts1, counts2), 
               samples = rbind(sampleinfo1, sampleinfo2))
print(paste0("** Dimensions of dge$counts: ", dim(dge$counts)))
print(paste0("** Dimensions of dge$samples: ", dim(dge$samples)))
print("** DGEList object: ")
print(dge)


# Variables for design matrix
sample_id <- factor(dge$samples$sample_name)
clusters_compare <- factor(dge$samples$cluster_compare)

# Define design matrix
design <- model.matrix(~ sample_id + clusters_compare)
colnames(design) <- make.names(colnames(design))
rownames(design) <- colnames(dge)

# filter the DGE object 
keep_genes <- filterByExpr(dge, design = design)
dge <- dge[keep_genes,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

# Estimate dispersion and fit the model
dge <- estimateDisp(dge, design)

# Perform DE testing
topRes <- dge %>% 
  glmQLFit(., design) %>%
  glmQLFTest(., coef = ncol(design)) %>%
  topTags(., n = Inf, p.value = 1) %>% 
  as.data.frame()

topRes <- topRes %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") 


# define output variables
this_comparison <- paste0("subcluster_6-2_vs_6-1")


# save the DE results
topRes %>%  
  arrange(PValue, -logFC) %>%
  write.csv(., file.path(opt$output, "gene_expression_associations",
                         paste0("pseudobulk_de_results_",
                                opt$output_prefix,
                                "_",
                                this_comparison, 
                                ".csv")), 
            row.names = FALSE)

# volcano plots
plot_volcano(topRes,
             cell_type = opt$output_prefix,
             comparison = this_comparison,
             p_value_cutoff = opt$p_val_cutoff,
             fold_change_cutoff = opt$fold_change_cutoff,
             output = file.path(opt$output, "gene_expression_associations")
)

# Plot lowest p-value box plots
plot_boxplots_top_genes(dge = dge, 
                        topRes = topRes, 
                        top_n = 50, 
                        output = file.path(opt$output, "gene_expression_associations"), 
                        filename = paste0("lowest_50_p_val_genes_boxplots_",
                                          opt$output_prefix,
                                          "_",
                                          this_comparison,
                                          ".pdf"))


# ranked genes for ORA
ranked_de_genes <- topRes %>%
  as.data.frame() %>%
  dplyr::select(gene, logFC, FDR) %>%
  mutate(rank = sign(logFC) * -log10(FDR)) %>%
  arrange(rank) %>%
  dplyr::rename(fold.change = logFC, p.value = FDR)


# ORA for ranked genes
run_ora(
  gene_list = ranked_de_genes,
  output = file.path(opt$output, "ORA"),
  output_prefix = paste0(opt$output_prefix,
                         "_",
                         this_comparison),
  fromType = "SYMBOL"
)


# GSEA for ranked genes
run_gsea(
  ranked_genes = ranked_de_genes,
  output = file.path(opt$output, "GSEA"),
  fromType = "SYMBOL",
  output_prefix = paste0(opt$output_prefix,
                         "_",
                         this_comparison)
)


# save the session info
writeLines(
  capture.output(sessionInfo()),
  file.path(opt$output, "sessionInfo.txt")
)


# END --------------------------------------------------------------------------