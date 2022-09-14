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
salmon_tpm_file <- snakemake@input[["salmon_tpm_file"]]

########
# MAIN #
########

mh_itWT_dds <- readRDS(mh_itWT_dds_file)

###############################
## iterative pairwise comp.s ##
###############################

## venom head
venom_head <- results(mh_itWT_dds, contrast=c("Tissue", "Venom", "Head"), alpha = 0.05, lfcThreshold = 1)
summary(venom_head)
## Order based of padj
ordered_venom_head <- venom_head[order(venom_head$padj),]
## Make data table and write to output
ordered_venom_head_table <- data.table(data.frame(ordered_venom_head), keep.rownames = TRUE)
head <- subset(ordered_venom_head_table, padj < 0.05)
fwrite(head, "output/03_deseq2/tissue_itWT/Venom/venom_vs_tissue_res/Venom_vs_head_res.csv")

## venom abdo
venom_abdo <- results(mh_itWT_dds, contrast=c("Tissue", "Venom", "Abdomen"), alpha = 0.05, lfcThreshold = 1)
summary(venom_abdo)
## Order based of padj
ordered_venom_abdo <- venom_abdo[order(venom_abdo$padj),]
## Make data table and write to output
ordered_venom_abdo_table <- data.table(data.frame(ordered_venom_abdo), keep.rownames = TRUE)
abdo <- subset(ordered_venom_abdo_table, padj < 0.05)
fwrite(abdo, "output/03_deseq2/tissue_itWT/Venom/venom_vs_tissue_res/Venom_vs_abdomen_res.csv")

## venom thorax
venom_thorax <- results(mh_itWT_dds, contrast=c("Tissue", "Venom", "Thorax"), alpha = 0.05, lfcThreshold = 1)
summary(venom_thorax)
## Order based of padj
ordered_venom_thorax <- venom_thorax[order(venom_thorax$padj),]
## Make data table and write to output
ordered_venom_thorax_table <- data.table(data.frame(ordered_venom_thorax), keep.rownames = TRUE)
thorax <- subset(ordered_venom_thorax_table, padj < 0.05)
fwrite(thorax, "output/03_deseq2/tissue_itWT/Venom/venom_vs_tissue_res/Venom_vs_thorax_res.csv")

## venom ovaries 
venom_ovaries <- results(mh_itWT_dds, contrast=c("Tissue", "Venom", "Ovaries"), alpha = 0.05, lfcThreshold = 1)
summary(venom_ovaries)
## Order based of padj
ordered_venom_ovaries <- venom_ovaries[order(venom_ovaries$padj),]
## Make data table and write to output
ordered_venom_ovaries_table <- data.table(data.frame(ordered_venom_ovaries), keep.rownames = TRUE)
ovaries <- subset(ordered_venom_ovaries_table, padj < 0.05)
fwrite(ordered_venom_ovaries_table, "output/03_deseq2/tissue_itWT/Venom/venom_vs_tissue_res/Venom_vs_ovaries_res.csv")

################################
## overlap for venom-specific ##
################################

## venn diagram
pdf(snakemake@output[["itWT_venn"]])
vd1 <- venn.diagram(x = list("Head vs Venom"=head$rn, "Thorax vs Venom"=thorax$rn, "Abomen vs Venom"=abdo$rn, "Ovaries vs Venom"=ovaries$rn), filename=NULL,
                    fill=c("#231151FF", "#5F187FFF", "#982D80FF", "#D3436EFF"), alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd1)
dev.off()

venom_specific_DEGs <- intersect(intersect(intersect(head$rn, thorax$rn), abdo$rn), ovaries$rn)
venom_v_ovary_ts <- subset(ovaries, rn %in% venom_specific_DEGs)

trinotate <- fread(trinotate_file, na.strings = ".")
venom_degs_annots <- merge(venom_v_ovary_ts, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(venom_degs_annots, snakemake@output[["tissue_annots"]])

##############
## DEG TPMs ##
##############

salmon_tpm <- fread(salmon_tpm_file)
summary_tpm <- data.table()
summary_tpm$rn <- salmon_tpm$rn
summary_tpm$venom_TPM <- rowMeans(salmon_tpm[,c("Mh_venom1", "Mh_venom2")], na.rm=TRUE)
summary_tpm$rest_TPM <- rowMeans(salmon_tpm[,!c("rn", "Mh_venom1", "Mh_venom2")], na.rm=TRUE)
setorder(summary_tpm, -venom_TPM)
summary_tpm$venom_rank <- seq.int(nrow(summary_tpm))
## merge with annotation table
annots_TPM <- merge(venom_degs_annots, summary_tpm, by="rn")
fwrite(annots_TPM, "output/03_deseq2/tissue_itWT/Venom/DEGs_annots_TPM.csv")

#############
## heatmap ##
#############

## vst transform
mh_vst <- varianceStabilizingTransformation(mh_itWT_dds, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)
##subset for DEGs
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% venom_specific_DEGs)
##turn first row back to row name
mh_vst_degs <- mh_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
## reorder tissues for plot - tissue of interest first
mh_vst_degs_plot <- mh_vst_degs[,c(13,14,4,5,6,16,10,11,12,1,2,3,15,7,8,9,17,18)]

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

# write log
sessionInfo()