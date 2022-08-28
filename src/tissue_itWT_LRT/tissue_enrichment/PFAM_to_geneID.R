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

library(data.table)
library(tidyverse)

###########
# GLOBALS #
###########

trinotate_file <- snakemake@input[["trinotate_file"]]

########
# MAIN #
########

trinotate_report <- fread(trinotate_file, na.strings=".")
##fairly sure we always used Pfam GO annots for this
pfam_annot_list<-data.table(trinotate_report[,unique(unlist(strsplit(Pfam, "`")))])
pfam_annot_table <- pfam_annot_list[,tstrsplit(V1, "^", fixed=TRUE)]
pfam_annot_table <- pfam_annot_table[,c(1,2,3)]
pfam_annot_table<-setnames(pfam_annot_table, old=c("V1", "V2", "V3"), new=c("Pfam_ID", "domain", "name"))
fwrite(pfam_annot_table, snakemake@output[["term_annot_table"]])

##function to extract GO terms from annotations in transcriptome (get all unique GO terms for each gene id) --> could look at other functional annot if I want to
EXTRACT_pfam_TERMS <- function(x, trinotate_report){
  my_terms<-trinotate_report[`#gene_id`==x,unique(unlist(strsplit(Pfam, "`")))]
  my_accessions<-unique(gsub("\\^.*", "", my_terms))
  my_accessions<-my_accessions[!is.na(my_accessions)]
  return(data.table(gene_id=x, accessions=my_accessions))
}

##extract GO terms for isoforms of genes
gene_ids <- trinotate_report[!is.na(Pfam), unique(`#gene_id`)]
pfam_term_list <- lapply(gene_ids, EXTRACT_pfam_TERMS, trinotate_report=trinotate_report)
pfam_term_table <- rbindlist(pfam_term_list)
##table of GO term to gene name
term_to_gene <- pfam_term_table[,c(2,1)]
term_to_name <- pfam_annot_table[,c(1,3)]
fwrite(distinct(term_to_gene), snakemake@output[["term_to_gene"]])
fwrite(distinct(term_to_name), snakemake@output[["term_to_name"]])

