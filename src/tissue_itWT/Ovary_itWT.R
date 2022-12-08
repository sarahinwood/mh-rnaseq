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
library(DESeq2)
library(VennDiagram)
library(pheatmap)
library(tidyverse)
library(viridis)

###########
# GLOBALS #
###########

mh_itWT_dds_file <- snakemake@input[["mh_itWT_dds_file"]]
trinotate_file <- snakemake@input[["trinotate_file"]]

########
# MAIN #
########

mh_itWT_dds <- readRDS(mh_itWT_dds_file)

###############################
## iterative pairwise comp.s ##
###############################

## ovaries head
ov_head <- results(mh_itWT_dds, contrast=c("Tissue", "Ovaries", "Head"), alpha = 0.05, lfcThreshold = 1)
summary(ov_head)
## Order based of padj
ordered_ov_head <- ov_head[order(ov_head$padj),]
## Make data table and write to output
ordered_ov_head_table <- data.table(data.frame(ordered_ov_head), keep.rownames = TRUE)
head <- subset(ordered_ov_head_table, padj < 0.05)
ordered_ov_head_table$tissue <- paste("head")

## ovaries abdo
ov_abdo <- results(mh_itWT_dds, contrast=c("Tissue", "Ovaries", "Abdomen"), alpha = 0.05, lfcThreshold = 1)
summary(ov_abdo)
## Order based of padj
ordered_ov_abdo <- ov_abdo[order(ov_abdo$padj),]
## Make data table and write to output
ordered_ov_abdo_table <- data.table(data.frame(ordered_ov_abdo), keep.rownames = TRUE)
abdo <- subset(ordered_ov_abdo_table, padj < 0.05)
ordered_ov_abdo_table$tissue <- paste("abdo")

## ovaries thorax
ov_thorax <- results(mh_itWT_dds, contrast=c("Tissue", "Ovaries", "Thorax"), alpha = 0.05, lfcThreshold = 1)
summary(ov_thorax)
## Order based of padj
ordered_ov_thorax <- ov_thorax[order(ov_thorax$padj),]
## Make data table and write to output
ordered_ov_thorax_table <- data.table(data.frame(ordered_ov_thorax), keep.rownames = TRUE)
thorax <- subset(ordered_ov_thorax_table, padj < 0.05)
ordered_ov_thorax_table$tissue <- paste("thorax")

## ovaries venom 
ov_venom <- results(mh_itWT_dds, contrast=c("Tissue", "Ovaries", "Venom"), alpha = 0.05, lfcThreshold = 1)
summary(ov_venom)
## Order based of padj
ordered_venom_ov <- ov_venom[order(ov_venom$padj),]
## Make data table and write to output
ordered_venom_ov_table <- data.table(data.frame(ordered_venom_ov), keep.rownames = TRUE)
venom <- subset(ordered_venom_ov_table, padj < 0.05)
ordered_venom_ov_table$tissue <- paste("venom")

tissue_res <- full_join(ordered_ov_head_table, full_join(ordered_ov_abdo_table, full_join(ordered_ov_thorax_table, ordered_venom_ov_table)))
fwrite(tissue_res, "output/03_deseq2/tissue_itWT/Ovary/Ovary_all_comps_results.csv")

################################
## overlap for ovary-specific ##
################################

## venn diagram
pdf(snakemake@output[["itWT_venn"]])
vd1 <- venn.diagram(x = list("Head"=head$rn, "Thorax"=thorax$rn, "Abomen"=abdo$rn, "Venom"=venom$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd1)
dev.off()

ovary_specific_DEGs <- intersect(intersect(intersect(head$rn, thorax$rn), abdo$rn), venom$rn)
ovary_v_venom_ts <- subset(venom, rn %in% ovary_specific_DEGs)

trinotate <- fread(trinotate_file, na.strings = ".")
ovary_degs_annots <- merge(ovary_v_venom_ts, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(ovary_degs_annots, snakemake@output[["tissue_annots"]])

#############
## heatmap ##
#############

## vst transform
mh_vst <- varianceStabilizingTransformation(mh_itWT_dds, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)
## subset for DEGs
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% ovary_specific_DEGs)
## turn first row back to row name
mh_vst_degs <- mh_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
## reorder tissues for plot - tissue of interest first
mh_vst_degs_plot <- mh_vst_degs[,c(7,8,9,17,18,4,5,6,16,10,11,12,1,2,3,15,13,14)]

## get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mh_itWT_dds)[,c("Tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
## for plot label
sample_to_tissue <- as.data.frame(colData(mh_itWT_dds)[,c("Tissue", "sample_name")])
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")

tissue_colours <- list(Tissue = c(Head='#2D1160FF', Thorax='#721F81FF', Abdomen='#B63679FF', Ovaries='#F1605DFF', Venom='#FEAF77FF'))

## plot
# not clustered by sample
pdf(snakemake@output[["heatmap"]])
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=FALSE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
dev.off()

# clustered by sample
pdf(snakemake@output[["clustered_heatmap"]])
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
dev.off()
# cutting off abdomen label from graph - open PDF in inkscape and resize page

# write log
sessionInfo()