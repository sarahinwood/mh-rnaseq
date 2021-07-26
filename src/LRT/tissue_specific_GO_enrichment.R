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
term_to_gene <- fread("output/deseq2/GO_annots/GO_to_geneID.csv")
term_to_name<- fread("output/deseq2/GO_annots/GO_to_GOname.csv")
go_annot_table <- fread("output/deseq2/GO_annots/GO_annots.csv")
go_annot_table <- go_annot_table[,c(1,2)]

results <- enricher(tissue_DEGs$rn, TERM2GENE=term_to_gene, TERM2NAME = term_to_name)
##seems legit enough?
results_dt <- data.table(as.data.frame(results), keep.rownames = FALSE)
##merge with annots table
plot_results_dt <- merge(results_dt, go_annot_table, by.x="ID", by.y="pathway")
plot_results_dt$total_DEG_GOs <- tstrsplit(plot_results_dt$GeneRatio, "/", keep=c(2))
plot_results_dt$total_DEG_GOs <- as.numeric(plot_results_dt$total_DEG_GOs)
plot_results_dt$GeneRatio <- (plot_results_dt$Count/plot_results_dt$total_DEG_GOs)
fwrite(plot_results_dt, snakemake@output[["enrichment_table"]])

## GeneRatio - number of DEGs with that GO annot/total DEGs with GO annots - DEG list
## BgRatio - total no. genes with that GO annot/total number of genes with GO annots - whole transcriptome
## Q value?

#####################
## all in one plot ##
#####################
###swap _ for space in pathway kind
plot_results_dt$pathway_kind <- gsub("_", " ", plot_results_dt$pathway_kind)
plot_results_dt$Description <- tstrsplit(plot_results_dt$Description, ",", keep=c(1))
##reorder - sorts by pathway_kind reverse alphabetically but can't figure out how to do any better
plot_results_dt$Description <- factor(plot_results_dt$Description, levels=plot_results_dt$Description[order(plot_results_dt$pathway_kind, plot_results_dt$GeneRatio, decreasing=TRUE)])

##plot
pdf(snakemake@output[['GO_plot']], height=3, width=5.875)
ggplot(plot_results_dt, aes(Description, GeneRatio)) +
  geom_col(aes(fill=pathway_kind))+
  labs(x="Gene ontology terms", y="GeneRatio",
       fill="GO domain", size="Leading\nedge size") +
  coord_flip() +
  scale_fill_viridis(discrete=TRUE)+
  theme_bw()
dev.off()

# write log
sessionInfo()
