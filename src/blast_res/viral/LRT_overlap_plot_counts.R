library(data.table)
library(DESeq2)
library(VennDiagram)

## LRT DEGs
mh_dds_lrt <- readRDS("output/03_deseq2/tissue_itWT_LRT/Venom/Venom_LRT_dds.rds")
res_group <- results(mh_dds_lrt, alpha = 0.05)
## Make data table and write to output
ordered_res_group_table <- data.table(data.frame(res_group), keep.rownames = TRUE)
lrt_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)

## Viral DEGs
fread("data/mh-transcriptome/")




pupa_dds <- readRDS("output/deseq2/stage_WT/mh_stage_WT.rds")
##prodigal viral
prodigal_viral_genes_blast <- fread("output/blast/viral_genes/viral_genes_best_hits.csv")
prodigal_viral_genes_blast$Trinity_ID <- tstrsplit(prodigal_viral_genes_blast$Trinity_ID,"_i", keep=c(1))
##crawford venom
crawford_venom_blast <- fread("output/blast/crawford_venom/crawford_best_hits.csv")
crawford_venom_blast$Trinity_ID <- tstrsplit(crawford_venom_blast$Trinity_ID,"_i", keep=c(1))

##overlap of LRT DEGs, venom, viral genes
vd1 <- venn.diagram(x = list("LRT DEGs"=LRT_sig_degs$rn, "Prodigal Viral"=prodigal_viral_genes_blast$Trinity_ID,
                               "Crawford venom"=crawford_venom_blast$Trinity_ID), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd1)

viral_degs <- intersect(LRT_sig_degs$rn, prodigal_viral_genes_blast$Trinity_ID)

###############################
## plot multiple gene counts ##
###############################

##get gene counts
counts_table <- data.table(counts(pupa_dds, normalized=TRUE), keep.rownames = TRUE)

## filter by list of interest - swap between crawford and prodigal to make both plots ##
annot_counts <- filter(counts_table, rn %in% prodigal_viral_genes_blast$Trinity_ID)

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
