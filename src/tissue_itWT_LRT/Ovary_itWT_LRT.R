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
mh_adult_dds_file <- snakemake@input[["mh_adult_dds_file"]]
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
head$tissue <- paste("head")

## ovaries abdo
ov_abdo <- results(mh_itWT_dds, contrast=c("Tissue", "Ovaries", "Abdomen"), alpha = 0.05, lfcThreshold = 1)
summary(ov_abdo)
## Order based of padj
ordered_ov_abdo <- ov_abdo[order(ov_abdo$padj),]
## Make data table and write to output
ordered_ov_abdo_table <- data.table(data.frame(ordered_ov_abdo), keep.rownames = TRUE)
abdo <- subset(ordered_ov_abdo_table, padj < 0.05)
abdo$tissue <- paste("abdo")

## ovaries thorax
ov_thorax <- results(mh_itWT_dds, contrast=c("Tissue", "Ovaries", "Thorax"), alpha = 0.05, lfcThreshold = 1)
summary(ov_thorax)
## Order based of padj
ordered_ov_thorax <- ov_thorax[order(ov_thorax$padj),]
## Make data table and write to output
ordered_ov_thorax_table <- data.table(data.frame(ordered_ov_thorax), keep.rownames = TRUE)
thorax <- subset(ordered_ov_thorax_table, padj < 0.05)
thorax$tissue <- paste("thorax")

## ovaries venom 
ov_venom <- results(mh_itWT_dds, contrast=c("Tissue", "Ovaries", "Venom"), alpha = 0.05, lfcThreshold = 1)
summary(ov_venom)
## Order based of padj
ordered_venom_ov <- ov_venom[order(ov_venom$padj),]
## Make data table and write to output
ordered_venom_ov_table <- data.table(data.frame(ordered_venom_ov), keep.rownames = TRUE)
venom <- subset(ordered_venom_ov_table, padj < 0.05)
venom$tissue <- paste("venom")

tissue_res <- full_join(head, full_join(abdo, full_join(thorax, venom)))

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

#####################
## LRT tissue test ##
#####################

mh_dds_lrt <- readRDS(mh_adult_dds_file)
## re-level tissue factor - ovaries are tissue of interest so put as first level for lfc calculation
mh_dds_lrt$Tissue <- factor(mh_dds_lrt$tissue, levels=c("Ovaries", "Head", "Thorax", "Abdomen", "Venom"))
mh_dds_lrt$Rep <- factor(mh_dds_lrt$rep)
mh_dds_lrt$Flowcell <- factor(mh_dds_lrt$flowcell)
## test
design(mh_dds_lrt) <- ~Flowcell+Rep+Tissue
mh_dds_lrt <- DESeq(mh_dds_lrt, test="LRT", reduced=~Flowcell+Rep)
## results
res_group <- results(mh_dds_lrt, alpha = 0.05)
summary(res_group)
## Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
## Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
lrt_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(lrt_sig_res_group_table, "output/03_deseq2/tissue_itWT_LRT/Ovary/sig_Ovary_lrt.csv")
saveRDS(mh_dds_lrt, "output/03_deseq2/tissue_itWT_LRT/Ovary/Ovary_LRT_dds.rds")

#########################
## overlap WT and LRT  ##
#########################

trinotate <- fread(trinotate_file, na.strings = ".")
LRT_ovary_degs <- subset(lrt_sig_res_group_table, rn %in% ovary_specific_DEGs)
LRT_ovary_degs_annots <- merge(LRT_ovary_degs, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(LRT_ovary_degs_annots, snakemake@output[["tissue_annots"]])

# unann_blast_res <- fread("output/03_deseq2/tissue_itWT_LRT/unann_degs/min_evalues.csv")
# LRT_ovary_degs_annots_blast <- merge(LRT_ovary_degs_annots, unann_blast_res, by="transcript_id", all.x=TRUE)

pdf(snakemake@output[["itWT_LRT_venn"]])
vd_ov <- venn.diagram(x = list("itWT Ovary DEGs"=ovary_specific_DEGs, "LRT"=lrt_sig_res_group_table$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_ov)
dev.off()

#############
## heatmap ##
#############

## vst transform
mh_vst <- varianceStabilizingTransformation(mh_dds_lrt, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)
## subset for DEGs
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% ovary_specific_DEGs)
## turn first row back to row name
mh_vst_degs <- mh_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
## reorder tissues for plot - tissue of interest first
mh_vst_degs_plot <- mh_vst_degs[,c(7,8,9,17,18,4,5,6,16,10,11,12,1,2,3,15,13,14)]

## get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
## for plot label
sample_to_tissue <- as.data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")])
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
