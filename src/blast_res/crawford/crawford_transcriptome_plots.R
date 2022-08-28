library(data.table)
library(tidyverse)
library(pheatmap)
library(viridis)
library(DESeq2)

############################
## heatmap venom proteins ##
############################

min_evalues_annots <- fread("output/04_blast/crawford_transcriptome/transcripts_best_blastn.csv")
##id to venom protein for plot
min_evalues_annots$label <- tstrsplit(min_evalues_annots$crawford_gene, "protein_", keep=c(2))
min_evalues_annots$label <- tstrsplit(min_evalues_annots$label, "_mRNA", keep=c(1))
min_evalues_annots$plot_label <- paste("Venom protein", min_evalues_annots$label, sep=" ")
id_venom <- min_evalues_annots[,c(13,34)]

mh_dds_lrt <- readRDS("output/03_deseq2/stage_WT/mh_stage_dds.rds")
##vst transform
mh_vst <- varianceStabilizingTransformation(mh_dds_lrt, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)

##subset for genes of interest
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% id_venom$`#gene_id`)
##merge with label
mh_vst_degs_label <- merge(mh_vst_degs, id_venom, by.x="rn", by.y="#gene_id.x")
##remove Trinity ID column
mh_vst_degs_label <- mh_vst_degs_label[,-c(1)]
setorder(mh_vst_degs_label, plot_label)
##turn first row back to row name
mh_vst_degs_label <- mh_vst_degs_label %>% remove_rownames %>% column_to_rownames(var="plot_label")
##reorder for plot
mh_vst_degs_plot <- mh_vst_degs_label[,c(16,17,20,8,7,9,21,4,5,6,19,13,14,15,1,2,3,18,10,11,12)]
mh_vst_degs_plot <- mh_vst_degs_plot[c(1,3,4,5,6,7,8,2),]

##get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mh_dds_lrt)[,c("tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
tissue_colours <- list(tissue = c(Head="#231151FF", Thorax="#5F187FFF", Abdomen="#982D80FF", Ovaries="#D3436EFF", Venom="#F8765CFF", Pupa="#FCFDBFFF"))
##plot
##not clustered by sample
pheatmap(mh_vst_degs_plot, cluster_rows=F, cluster_cols=F, show_rownames=TRUE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))

# Modify ordering of the clusters using clustering callback option - pheatmap manual
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
##clustered by gene
pheatmap(mh_vst_degs_plot, cluster_rows=T, cluster_cols=T, clustering_callback=callback,
         show_rownames=TRUE, cutree_cols = 4,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))

## what about without ovary and other abdo sample?
mh_dds_lrt_subset <- mh_dds_lrt[, mh_dds_lrt$Flowcell=="HM723BCX3"]
mh_vst_s <- varianceStabilizingTransformation(mh_dds_lrt_subset, blind=TRUE)
mh_vst_assay_dt_s <- data.table(assay(mh_vst_s), keep.rownames=TRUE)
##subset for genes of interest
mh_vst_degs_s <- subset(mh_vst_assay_dt_s, rn %in% id_venom$`#gene_id`)
##merge with label
mh_vst_degs_label_s <- merge(mh_vst_degs_s, id_venom, by.x="rn", by.y="#gene_id.x")
##remove Trinity ID column
mh_vst_degs_label_s <- mh_vst_degs_label_s[,-c(1)]
setorder(mh_vst_degs_label_s, plot_label)
##turn first row back to row name
mh_vst_degs_label_s <- mh_vst_degs_label_s %>% remove_rownames %>% column_to_rownames(var="plot_label")
##reorder for plot
#mh_vst_degs_plot_s <- mh_vst_degs_label[,c(16,17,20,8,7,9,21,4,5,6,19,13,14,15,1,2,3,18,10,11,12)]
#mh_vst_degs_plot <- mh_vst_degs_plot[c(1,3,4,5,6,7,8,2),]

##get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mh_dds_lrt_subset)[,c("tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
tissue_colours <- list(tissue = c(Head="#231151FF", Thorax="#5F187FFF", Abdomen="#982D80FF", Ovaries="#D3436EFF", Venom="#F8765CFF", Pupa="#FCFDBFFF"))

# Modify ordering of the clusters using clustering callback option - pheatmap manual
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
##clustered by gene
pheatmap(mh_vst_degs_label_s, cluster_rows=T, cluster_cols=T, clustering_callback=callback,
         show_rownames=TRUE, cutree_cols = 3,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))

  
############################
## plot venom gene counts ##
############################
##plot venom genes not DE
mh_dds_lrt <- readRDS("output/deseq2/stage_WT/mh_stage_WT.rds")
##for best hits##
counts_table <- data.table(counts(mh_dds_lrt, normalized=TRUE), keep.rownames = TRUE)
best_annot_counts <- filter(counts_table, rn %in% min_evalues_annots$`#gene_id`)
##melt for plotting
plot_annots_counts <- best_annot_counts %>% gather(colnames(best_annot_counts)[2:22], key="sample_name", value="normalized_counts")
##sample group information
sample_table <- fread("data/sample_table.csv")
sample_table$group <- paste(sample_table$tissue)
name_vs_group <- sample_table[,c(1,2)]
plotting_counts <- inner_join(plot_annots_counts, name_vs_group)
group_order <- c("Head", "Thorax", "Abdomen", "Ovaries", "Venom", "Pupa")
plotting_counts$group <- factor(plotting_counts$tissue, levels=group_order)
##merge with crawford gene annot
min_evalues_annots$label <- tstrsplit(min_evalues_annots$Crawford_seq, "protein_", keep=c(2))
min_evalues_annots$label <- tstrsplit(min_evalues_annots$label, "_mRNA", keep=c(1))
min_evalues_annots$plot_label <- paste("Venom protein", min_evalues_annots$label, sep=" ")
vp_order <- c("Venom protein 1", "Venom protein 2", "Venom protein 3","Venom protein 4","Venom protein 5", "Venom protein 6", "Venom protein 8", "Venom protein 10")
min_evalues_annots$plot_label <- factor(min_evalues_annots$plot_label, levels=vp_order)
gene_to_crawford <- min_evalues_annots[,c(14,31)]
plotting_counts_crawford_hits <- merge(plotting_counts, gene_to_crawford, by.x="rn", by.y="#gene_id")
##plot all annot DEGs using ggplot2
ggplot(plotting_counts_crawford_hits) +
  geom_point(aes(x = group, y = normalized_counts, colour=group)) +
  labs(colour="Tissue", y="Normalized counts", x="")+
  scale_colour_viridis(discrete=TRUE)+
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~plot_label, scales="free")

##for all jhits
##get gene counts
annot_counts <- filter(counts_table, rn %in% crawford_blast_res$`#gene_id`)
##melt for plotting
plot_annots_counts <- annot_counts %>% gather(colnames(annot_counts)[2:22], key="sample_name", value="normalized_counts")
##sample group information
sample_table <- fread("data/sample_table.csv")
sample_table$group <- paste(sample_table$tissue)
name_vs_group <- sample_table[,c(1,2)]
plotting_counts <- inner_join(plot_annots_counts, name_vs_group)
group_order <- c("Head", "Thorax", "Abdomen", "Ovaries", "Venom", "Pupa")
plotting_counts$group <- factor(plotting_counts$tissue, levels=group_order)
##merge with crawford gene annot
crawford_blast_res$label <- tstrsplit(crawford_blast_res$Crawford_seq, "hyperodae_", keep=c(2))
crawford_blast_res$label <- tstrsplit(crawford_blast_res$label, "_mRNA", keep=c(1))
gene_to_crawford <- crawford_blast_res[,c(14,15,16)]
plotting_counts_crawford_hits <- merge(plotting_counts, gene_to_crawford, by.x="rn", by.y="#gene_id")
plotting_counts_crawford_hits$plot_label <- paste(plotting_counts_crawford_hits$short_gene_id, plotting_counts_crawford_hits$label, sep="_")
##plot all annot DEGs using ggplot2
ggplot(plotting_counts_crawford_hits) +
  geom_point(aes(x = group, y = normalized_counts, colour=group)) +
  labs(colour="Tissue", y="Normalized counts", x="")+
  scale_colour_viridis(discrete=TRUE)+
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~plot_label, scales="free")

