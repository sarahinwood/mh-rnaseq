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
pfam_annot_table_file <- snakemake@input[["term_annot_table_file"]]

########
# MAIN #
########

tissue_DEGs <- fread(tissue_DEG_file)

##GO annot tables
pfam_to_gene <- fread(term_to_gene_file)
pfam_to_name<- distinct(fread(term_to_name_file))
pfamannot_table <- fread(pfam_annot_table_file)
pfamannot_table <- pfamannot_table[,c(1,2)]
pfamannot_table <- distinct(pfamannot_table)

results <- enricher(tissue_DEGs$rn, TERM2GENE=pfam_to_gene, TERM2NAME = pfam_to_name)
##seems legit enough?
results_dt <- data.table(as.data.frame(results), keep.rownames = FALSE)
##merge with annots table
plot_results_dt <- merge(results_dt, pfamannot_table, by.x="ID", by.y="Pfam_ID", all.x=TRUE)
plot_results_dt$total_DEG_PFAMs <- tstrsplit(plot_results_dt$GeneRatio, "/", keep=c(2))
plot_results_dt$total_DEG_PFAMs <- as.numeric(plot_results_dt$total_DEG_PFAMs)
plot_results_dt$GeneRatio <- (plot_results_dt$Count/plot_results_dt$total_DEG_PFAMs)

fwrite(plot_results_dt, snakemake@output[["enrichment_table"]])

# write log
sessionInfo()
