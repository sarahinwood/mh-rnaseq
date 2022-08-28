library(data.table)
library(ggplot2)
library(viridis)
library(tidyverse)
library(DESeq2)

pupa_dds <- readRDS("output/deseq2/stage_WT/mh_stage_WT.rds")
LRT_sig_degs <- fread("output/deseq2/tissue_LRT/sig_degs.csv")

LbFV_DaFV_genes <- fread("data/mh-transcriptome/output/recip_blast/nr_blastx/LbFV_DaFV_annots.csv")
LbFV_DaFV_genes$`#gene_id` <- tstrsplit(LbFV_DaFV_genes$transcript_id, "_i", keep=c(1))
##sig DEG LbFV/DaFV
intersect <- intersect(LbFV_DaFV_genes$`#gene_id`, LRT_sig_degs$rn)

###############################
## plot multiple gene counts ##
###############################

##get gene counts
counts_table <- data.table(counts(pupa_dds, normalized=TRUE), keep.rownames = TRUE)

## filter by list of interest - swap between crawford and prodigal to make both plots ##
annot_counts <- filter(counts_table, rn %in% LbFV_DaFV_genes$`#gene_id`)

##melt for plotting
plot_annots_counts <- annot_counts %>% gather(colnames(annot_counts)[2:18], key="sample_name", value="normalized_counts")
##sample group information
sample_table <- fread("data/sample_table.csv")
name_vs_group <- sample_table[,c(1,2)]
plotting_counts <- inner_join(plot_annots_counts, name_vs_group)
tissue_order <- c("Head", "Thorax", "Abdomen", "Ovaries", "Venom", "Pupa")
plotting_counts$tissue <- factor(plotting_counts$tissue, levels=tissue_order)
##add alphabetical label to each plot
plotting_counts$gene_label <- paste(plotting_counts$rn)
##plot all annot DEGs using ggplot2
ggplot(plotting_counts) +
  geom_point(aes(x = tissue, y = normalized_counts, colour=tissue)) +
  labs(colour="Tissue", y="Normalized counts", x="")+
  scale_colour_viridis(discrete=TRUE)+
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~gene_label, scales="free")
