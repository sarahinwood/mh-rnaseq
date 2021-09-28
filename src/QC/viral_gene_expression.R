library(data.table)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

mh_dds <- readRDS("output/deseq2/tissue_LRT/mh_tissue_LRT.rds")
mh_viral <- fread("data/mh-transcriptome/output/trinotate/prokaryotic/viral_genes_plot.csv")
gene_family <- mh_viral[,c(1,10)]
high_rep_family <- list("Siphoviridae", "Retroviridae", "Rhabdoviridae", "Baculoviridae")
gene_family_oi <- subset(gene_family, family %in% high_rep_family)

###############################
## plot multiple gene counts ##
###############################

##get gene counts
counts_table <- data.table(counts(mh_dds, normalized=TRUE), keep.rownames = TRUE)

## filter by list of interest - swap between crawford and prodigal to make both plots ##
annot_counts <- filter(counts_table, rn %in% gene_family_oi$`#gene_id`)

##melt for plotting
plot_annots_counts <- annot_counts %>% gather(colnames(annot_counts)[2:15], key="sample_name", value="normalized_counts")
##merge with list of phyla
plot_annots_counts <- merge(plot_annots_counts, gene_family, by.x="rn", by.y="#gene_id")

##sample group information
sample_table <- fread("data/sample_table.csv")
name_vs_group <- sample_table[,c(1,2)]
plotting_counts <- inner_join(plot_annots_counts, name_vs_group)
tissue_order <- c("head", "thorax", "abdomen", "ovaries", "venom")
plotting_counts$tissue <- factor(plotting_counts$tissue, levels=tissue_order)

##plot all annot DEGs using ggplot2
ggplot(plotting_counts) +
  geom_point(aes(x = tissue, y = normalized_counts, colour=tissue)) +
  labs(colour="Tissue", y="Normalized counts", x="")+
  scale_colour_viridis(discrete=TRUE)+
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~family, scales="free")
