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

###########
# GLOBALS #
###########

trinotate_file <- snakemake@input[["trinotate_file"]]

########
# MAIN #
########

trinotate_report <- fread(trinotate_file, na.strings=".")
gene_ids <- trinotate_report[!is.na(gene_ontology_Pfam), unique(`#gene_id`)]
##fairly sure we always used Pfam GO annots for this
go_annot_list<-data.table(trinotate_report[,unique(unlist(strsplit(gene_ontology_Pfam, "`")))])
go_annot_table <- go_annot_list[,tstrsplit(V1, "^", fixed=TRUE)]
go_annot_table<-setnames(go_annot_table, old=c("V1", "V2", "V3"), new=c("pathway", "pathway_kind", "pathway_name"))
fwrite(go_annot_table, snakemake@output[["term_annot_table"]])

##function to extract GO terms from annotations in transcriptome (get all unique GO terms for each gene id) --> could look at other functional annot if I want to
EXTRACT_GO_TERMS <- function(x, trinotate_report){
  my_terms<-trinotate_report[`#gene_id`==x,unique(unlist(strsplit(gene_ontology_Pfam, "`")))]
  my_accessions<-unique(gsub("\\^.*", "", my_terms))
  my_accessions<-my_accessions[!is.na(my_accessions)]
  return(data.table(gene_id=x, accessions=my_accessions))
}

##extract GO terms for isoforms of genes
go_term_list <- lapply(gene_ids, EXTRACT_GO_TERMS, trinotate_report=trinotate_report)
go_term_table <- rbindlist(go_term_list)
##table of GO term to gene name
term_to_gene <- go_term_table[,c(2,1)]
term_to_name <- go_annot_table[,c(1,3)]
fwrite(term_to_gene, snakemake@output[["term_to_gene"]])
fwrite(term_to_name, snakemake@output[["term_to_name"]])

