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
library(tidyverse)
library(pheatmap)
library(viridis)
library(dendsort)

#############
## GLOBALS ##
#############

ovary_itWT_file <- snakemake@input[["ovary_itWT_file"]]
dds_file <- snakemake@input[["dds_file"]]
meiosis_sp_genes_file <- snakemake@input[["meiosis_sp_genes_file"]]
blastx_res <- snakemake@input[["blastx_res"]]

########
# MAIN #
########
meiosis_sp_blast <- fread(blastx_res)
meiosis_sp_genes <- fread(meiosis_sp_genes_file)

# remove short redundant transcripts for REC8 and DMC1
meiosis_sp_genes_fil <- subset(meiosis_sp_genes, !grepl("TRINITY_DN85877_c0_g1_i1|TRINITY_DN99454_c0_g1_i1|TRINITY_DN7176_c0_g1_i1", transcript_id))

lrt_degs <- fread(ovary_itWT_file)
meiosis_lrt <- subset(lrt_degs, lrt_degs$rn %in% meiosis_sp_genes_fil$`#gene_id`)
fwrite(meiosis_lrt, "output/03_deseq2/tissue_itWT_LRT/meiosis_sp/msp_LRT_res.csv")
##label for DE in LRT
meiosis_sp_genes_fil$DE <- ifelse(meiosis_sp_genes_fil$`#gene_id` %in% lrt_degs$rn, "*", "")
meiosis_sp_genes_fil$Gene <- paste(meiosis_sp_genes_fil$meiosis_gene, meiosis_sp_genes_fil$DE, sep="|")

gene_id_to_gene <- meiosis_sp_genes_fil[,c(2,23)]
setorder(gene_id_to_gene, Gene)
gene_id_to_gene$gene_id_gene <- paste(gene_id_to_gene$Gene, gene_id_to_gene$`#gene_id`, sep="_")
label <- gene_id_to_gene[,c(1,3)]

#############
## heatmap ##
#############

mh_dds_lrt <- readRDS(dds_file)
##vst transform
mh_vst <- varianceStabilizingTransformation(mh_dds_lrt, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)

##subset for genes of interest
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% gene_id_to_gene$`#gene_id`)
##merge with label
mh_vst_degs_label <- merge(mh_vst_degs, label, by.x="rn", by.y="#gene_id")
mh_vst_degs_label <- mh_vst_degs_label[,-c(1)]
setorder(mh_vst_degs_label, gene_id_gene)
##turn first row back to row name
mh_vst_degs_label <- mh_vst_degs_label %>% remove_rownames %>% column_to_rownames(var="gene_id_gene")

##get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
##for plot label
sample_to_tissue <- as.data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")])
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")

tissue_colours <- list(Tissue = c(Head="#231151FF", Thorax="#5F187FFF", Abdomen="#982D80FF", Ovaries="#D3436EFF", Venom="#F8765CFF"))
##plot ordering
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
##clustered by sample
pdf(snakemake@output[["heatmap"]], width=10, height=4)
pheatmap(mh_vst_degs_label, cluster_rows=T, cluster_cols=T, clustering_callback=callback, show_rownames=TRUE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50), cutree_cols = 2)
dev.off()

#open PDF in inkscape and change gene labels to remove trinity IDs, bold sig DEGs
##dendsort could be used to sort dendrogram? but may be not sorted by longest branch ue to sorting off expression patterns
