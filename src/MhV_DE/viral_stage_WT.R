library(data.table)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(viridis)
library(pheatmap)

viral_genes_blast <- fread("output/blast/viral_genes/viral_genes_best_hits.csv")
viral_genes_blast$Trinity_gene_ID <- tstrsplit(viral_genes_blast$Trinity_ID, "_i", keep=c(1))
viral_genes <- unique(viral_genes_blast$Trinity_gene_ID)

mh_dds <- readRDS("output/deseq2/mh_dds.rds")
##estimate size factors on whole table - differences in library size between samples
mh_dds <- estimateSizeFactors(mh_dds)

##keep only viral genes
mh_dds_viral <- mh_dds[viral_genes,]
##re-level factor
mh_dds_viral$Stage <- factor(mh_dds_viral$stage, levels=c("pupa", "adult"))
mh_dds_viral$Rep <- factor(mh_dds_viral$rep)

design(mh_dds_viral) <- ~Rep+Stage
mh_dds_viral <- DESeq(mh_dds_viral)
saveRDS(mh_dds_viral, "output/deseq2/MhV_LRT/mh_stage_LRT.rds")

res_group <- results(mh_dds_viral, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/deseq2/MhV_LRT/stage_sig_degs.csv")
fwrite(list(unique(sig_annots$rn)), "output/deseq2/MhV_LRT/stage_degs.txt")

trinotate <- fread("data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings = ".")
sig_trinotate <- merge(ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id")

#############
## heatmap ##
#############

mh_dds_viral <- readRDS("output/deseq2/MhV_LRT/mh_stage_LRT.rds")
ordered_sig_res_group_table <- fread("output/deseq2/MhV_LRT/stage_sig_degs.csv")

##vst transform
mh_vst <- varianceStabilizingTransformation(mh_dds_viral, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)
##subset for DEGs
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% ordered_sig_res_group_table$rn)
##turn first row back to row name
mh_vst_degs <- mh_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
##reorder for plot
mh_vst_degs_plot <- mh_vst_degs[,c(10,11,12,4,5,6,19,13,14,15,1,2,3,18,7,8,9,20,21,16,17)]

##get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mh_dds_viral)[,c("tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
##for plot label
sample_to_tissue <- as.data.frame(colData(mh_dds_viral)[,c("tissue", "sample_name")])
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")

tissue_colours <- list(tissue = c(Pupa="#440154FF", Head="#414487FF", Thorax="#2A788EFF", Abdomen="#22A884FF", Ovaries="#7AD151FF", Venom="#FDE725FF"))

##clustered by sample
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=FALSE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
