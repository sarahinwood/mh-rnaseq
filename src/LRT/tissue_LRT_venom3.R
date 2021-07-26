library(tximport)
library(data.table)
library(DESeq2)
library(VennDiagram)

mh_dds <- readRDS("output/deseq2/mh_dds_all.rds")

##select only adult tissue samples
mh_dds <- mh_dds[,mh_dds$stage=="Adult"]
##re-level tissue factor
mh_dds$Tissue <- factor(mh_dds$tissue, levels=c("Head", "Thorax", "Abdomen", "Ovaries", "Venom"))
mh_dds$Rep <- factor(mh_dds$rep)

design(mh_dds) <- ~Rep+Tissue
mh_dds <- DESeq(mh_dds, test="LRT", reduced=~Rep)
saveRDS(mh_dds, "output/deseq2/tissue_LRT_venom3/mh_tissue_LRT.rds")

res_group <- results(mh_dds, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/deseq2/tissue_LRT_venom3/sig_degs.csv")

plotCounts(mh_dds, "TRINITY_DN5452_c0_g1", intgroup=("tissue"))

##prodigal viral
prodigal_viral_genes_blast <- fread("output/blast/viral_genes/viral_genes_best_hits.csv")
prodigal_viral_genes_blast$Trinity_ID <- tstrsplit(prodigal_viral_genes_blast$Trinity_ID,"_i", keep=c(1))
##crawford venom - does 3rd sample show venom expression as it should or is it a total loss?
crawford_venom_blast <- fread("output/blast/crawford_venom/crawford_best_hits.csv")
crawford_venom_blast$Trinity_ID <- tstrsplit(crawford_venom_blast$Trinity_ID,"_i", keep=c(1))

###############################
## plot multiple gene counts ##
###############################

##get gene counts
counts_table <- data.table(counts(mh_dds, normalized=TRUE), keep.rownames = TRUE)

## filter by list of interest - swap between crawford and prodigal to make both plots ##
annot_counts <- filter(counts_table, rn %in% prodigal_viral_genes_blast$Trinity_ID)

##melt for plotting
plot_annots_counts <- annot_counts %>% gather(colnames(annot_counts)[2:16], key="sample_name", value="normalized_counts")
##sample group information
sample_table <- fread("data/sample_table.csv")
name_vs_group <- sample_table[,c(1,2)]
plotting_counts <- inner_join(plot_annots_counts, name_vs_group)
tissue_order <- c("Head", "Thorax", "Abdomen", "Ovaries", "Venom")
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

##venom genes expressed but much lower
##viral genes not expressed at all

##cooks with venom3?