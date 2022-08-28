#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

#############
# LIBRARIES #
#############

library(tidyverse)
library(data.table)
library(clusterProfiler)
library(viridis)

###########
# GLOBALS #
###########

##DEG_list             
tissue_DEG_file <- snakemake@input[["tissue_DEG_file"]]
term_to_gene_file <- snakemake@input[["term_to_gene_file"]]
term_to_name_file <- snakemake@input[["term_to_name_file"]]
go_annot_table_file <- snakemake@input[["term_annot_table_file"]]

########
# MAIN #
########

tissue_DEGs <- fread(tissue_DEG_file)

##GO annot tables
term_to_gene <- fread(term_to_gene_file)
term_to_name<- fread(term_to_name_file)
go_annot_table <- fread(go_annot_table_file)
go_annot_table <- go_annot_table[,c(1,2)]

results <- enricher(tissue_DEGs$rn, TERM2GENE=term_to_gene, TERM2NAME = term_to_name)
##seems legit enough?
results_dt <- data.table(as.data.frame(results), keep.rownames = FALSE)
##merge with annots table
plot_results_dt <- merge(results_dt, go_annot_table, by.x="ID", by.y="pathway")
plot_results_dt$total_DEG_GOs <- tstrsplit(plot_results_dt$GeneRatio, "/", keep=c(2))
plot_results_dt$total_DEG_GOs <- as.numeric(plot_results_dt$total_DEG_GOs)
plot_results_dt$GeneRatio <- (plot_results_dt$Count/plot_results_dt$total_DEG_GOs)

plot_results_dt$pathway_kind <- gsub("_", " ", plot_results_dt$pathway_kind)
plot_results_dt$Description <- tstrsplit(plot_results_dt$Description, ",", keep=c(1))

fwrite(plot_results_dt, snakemake@output[["enrichment_table"]])

## GeneRatio - number of DEGs with that GO annot/total DEGs with GO annots - DEG list
## BgRatio - total no. genes with that GO annot/total number of genes with GO annots - whole transcriptome
## Q value?

# write log
sessionInfo()
