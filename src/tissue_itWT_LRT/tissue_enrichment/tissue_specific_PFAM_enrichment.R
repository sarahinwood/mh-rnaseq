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
tissue_DEG_file <- snakemake@input[["tissue_degs_file"]]

########
# MAIN #
########

tissue_DEGs <- fread(tissue_DEG_file)

##GO annot tables
pfam_to_gene <- fread("output/deseq2/pfam_annots/pfam_to_geneID.csv")
pfam_to_name<- distinct(fread("output/deseq2/pfam_annots/pfam_to_domain_name.csv"))
pfamannot_table <- fread("output/deseq2/pfam_annots/pfam_annots.csv")
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

## GeneRatio - number of DEGs with that GO annot/total DEGs with GO annots - DEG list
## BgRatio - total no. genes with that GO annot/total number of genes with GO annots - whole transcriptome
## Q value?

#####################
## all in one plot ##
#####################
plot_results_dt$Description <- factor(plot_results_dt$Description, levels=(plot_results_dt$Description[order(plot_results_dt$GeneRatio, decreasing=TRUE)]))

##plot
pdf(snakemake@output[['PFAM_plot']], height=3, width=5.875)
ggplot(plot_results_dt, aes(x=Description, y=GeneRatio)) +
  geom_col(aes(fill="#440154FF"))+
  labs(x="Pfam domain", y="GeneRatio", size="Leading\nedge size") +
  coord_flip() +
  scale_fill_viridis(discrete=TRUE)+
  theme_bw()+
  theme(legend.position = "none")
dev.off()

# write log
sessionInfo()
