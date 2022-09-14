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

## abdo head
abdo_head <- results(mh_itWT_dds, contrast=c("Tissue", "Abdomen", "Head"), alpha = 0.05, lfcThreshold = 1)
summary(abdo_head)
## Order based of padj
ordered_abdo_head <- abdo_head[order(abdo_head$padj),]
## Make data table and write to output
ordered_abdo_head_table <- data.table(data.frame(ordered_abdo_head), keep.rownames = TRUE)
head <- subset(ordered_abdo_head_table, padj < 0.05)

## abdo thorax
abdo_thorax <- results(mh_itWT_dds, contrast=c("Tissue", "Abdomen", "Thorax"), alpha = 0.05, lfcThreshold = 1)
summary(abdo_thorax)
## Order based of padj
ordered_abdo_thorax <- abdo_thorax[order(abdo_thorax$padj),]
## Make data table and write to output
ordered_abdo_thorax_table <- data.table(data.frame(ordered_abdo_thorax), keep.rownames = TRUE)
thorax <- subset(ordered_abdo_thorax_table, padj < 0.05)

## abdo thorax
abdo_venom <- results(mh_itWT_dds, contrast=c("Tissue", "Abdomen", "Venom"), alpha = 0.05, lfcThreshold = 1)
summary(abdo_venom)
## Order based of padj
ordered_abdo_venom <- abdo_venom[order(abdo_venom$padj),]
## Make data table and write to output
ordered_abdo_venom_table <- data.table(data.frame(ordered_abdo_venom), keep.rownames = TRUE)
venom <- subset(ordered_abdo_venom_table, padj < 0.05)

## abdo ovaries 
abdo_ovaries <- results(mh_itWT_dds, contrast=c("Tissue", "Abdomen", "Ovaries"), alpha = 0.05, lfcThreshold = 1)
summary(abdo_ovaries)
## Order based of padj
ordered_abdo_ovaries <- abdo_ovaries[order(abdo_ovaries$padj),]
## Make data table and write to output
ordered_abdo_ovaries_table <- data.table(data.frame(ordered_abdo_ovaries), keep.rownames = TRUE)
ovaries <- subset(ordered_abdo_ovaries_table, padj < 0.05)

##################################
## overlap for abdomen-specific ##
##################################

## venn diagram - abdomen-specific DEGs
pdf(snakemake@output[["itWT_venn"]])
vd1 <- venn.diagram(x = list("Head"=head$rn, "Venom"=venom$rn, "Thorax"=thorax$rn, "Ovaries"=ovaries$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd1)
dev.off()

abdo_specific_DEGs <- intersect(intersect(intersect(head$rn, venom$rn), thorax$rn), ovaries$rn)
abdo_v_ovary_ts <- subset(ovaries, rn %in% abdo_specific_DEGs)

trinotate <- fread(trinotate_file, na.strings = ".")
abdo_degs_annots <- merge(abdo_v_ovary_ts, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(abdo_degs_annots, snakemake@output[["tissue_annots"]])

#############
## heatmap ##
#############

## vst transform
mh_vst <- varianceStabilizingTransformation(mh_itWT_dds, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)
## subset for DEGs
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% abdo_specific_DEGs)
## turn first row back to row name
mh_vst_degs <- mh_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
## reorder tissues for plot - tissue of interest first
mh_vst_degs_plot <- mh_vst_degs[,c(1,2,3,15,4,5,6,16,10,11,12,7,8,9,17,18,13,14)]

## get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mh_itWT_dds)[,c("Tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
## for plot label
sample_to_tissue <- as.data.frame(colData(mh_itWT_dds)[,c("Tissue", "sample_name")])
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")

# sample colour labels
tissue_colours <- list(Tissue = c(Head='#2D1160FF', Thorax='#721F81FF', Abdomen='#B63679FF', Ovaries='#F1605DFF', Venom='#FEAF77FF'))
## plot
## not clustered by samples
pdf(snakemake@output[["heatmap"]])
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=FALSE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
dev.off()

## clustered by samples
pdf(snakemake@output[["clustered_heatmap"]])
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
dev.off()

# write log
sessionInfo()