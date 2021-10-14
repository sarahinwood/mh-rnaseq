library(data.table)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

mh_dds <- readRDS("output/deseq2/tissue_itWT_LRT/mh_itWT.rds")
mh_bacterial <- fread("data/mh-transcriptome/output/trinotate/prokaryotic/bacterial_genes.csv")
gene_phyla <- mh_bacterial[,c(1,10)]
high_rep_phyla <- list("Actinobacteria", "Firmicutes", "Proteobacteria", "Bacteroidetes")
gene_phyla_oi <- subset(gene_phyla, phyla %in% high_rep_phyla)

###############################
## plot multiple gene counts ##
###############################

##get gene counts
counts_table <- data.table(counts(mh_dds, normalized=TRUE), keep.rownames = TRUE)

## filter by list of interest - swap between crawford and prodigal to make both plots ##
annot_counts <- filter(counts_table, rn %in% gene_phyla_oi$`#gene_id`)

##melt for plotting
plot_annots_counts <- annot_counts %>% gather(colnames(annot_counts)[2:19], key="sample_name", value="normalized_counts")
##merge with list of phyla
plot_annots_counts <- merge(plot_annots_counts, gene_phyla_oi, by.x="rn", by.y="#gene_id")

##sample group information
sample_table <- fread("data/sample_table.csv")
name_vs_group <- sample_table[,c(1,2)]
plotting_counts <- inner_join(plot_annots_counts, name_vs_group)
tissue_order <- c("Head", "Thorax", "Abdomen", "Ovaries", "Venom")
plotting_counts$tissue <- factor(plotting_counts$tissue, levels=tissue_order)

##LRT tissue analysis
tissue_LRT <- fread("output/deseq2/tissue_LRT/sig_degs.csv")
plotting_counts_LRT <- subset(plotting_counts, rn %in% tissue_LRT$rn)

##plot all annot DEGs using ggplot2
ggplot(plotting_counts_LRT) +
  geom_point(aes(x = tissue, y = normalized_counts, colour=tissue)) +
  labs(colour="Tissue", y="Normalized counts", x="")+
  scale_colour_viridis(discrete=TRUE)+
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~phyla, scales="free", nrow = 1)

trinotate <- fread("data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv")
deg_annots <- subset(trinotate, `#gene_id` %in% plotting_counts_LRT$rn)
TPMs <- fread("output/deseq2/salmon_TPM_means.csv")
annots_tpms <- merge(deg_annots, TPMs, by.x="#gene_id", by.y="rn", all.x=TRUE)
fwrite(annots_tpms, "output/prokaryotic_phyla_counts/bacterial_deg_annots.csv")

