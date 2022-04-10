library(data.table)
library(DESeq2)
library(VennDiagram)
library(viridis)
library(pheatmap)
library(tidyverse)

mh_itWT_dds <- readRDS("output/deseq2/tissue_itWT_LRT/mh_itWT.rds")

###############################
## iterative pairwise comp.s ##
###############################

##thorax head
thorax_head <- results(mh_itWT_dds, contrast=c("Tissue", "Thorax", "Head"), alpha = 0.05, lfcThreshold = 1)
summary(thorax_head)
##Order based of padj
ordered_thorax_head <- thorax_head[order(thorax_head$padj),]
##Make data table and write to output
ordered_thorax_head_table <- data.table(data.frame(ordered_thorax_head), keep.rownames = TRUE)
head <- subset(ordered_thorax_head_table, padj < 0.05)

##thorax abdo
thorax_abdo <- results(mh_itWT_dds, contrast=c("Tissue", "Thorax", "Abdomen"), alpha = 0.05, lfcThreshold = 1)
summary(thorax_abdo)
##Order based of padj
ordered_thorax_abdo <- thorax_abdo[order(thorax_abdo$padj),]
##Make data table and write to output
ordered_thorax_abdo_table <- data.table(data.frame(ordered_thorax_abdo), keep.rownames = TRUE)
abdo <- subset(ordered_thorax_abdo_table, padj < 0.05)

##thorax thorax
thorax_venom <- results(mh_itWT_dds, contrast=c("Tissue", "Thorax", "Venom"), alpha = 0.05, lfcThreshold = 1)
summary(thorax_venom)
##Order based of padj
ordered_thorax_venom <- thorax_venom[order(thorax_venom$padj),]
##Make data table and write to output
ordered_thorax_venom_table <- data.table(data.frame(ordered_thorax_venom), keep.rownames = TRUE)
venom <- subset(ordered_thorax_venom_table, padj < 0.05)

##thorax ovaries 
thorax_ovaries <- results(mh_itWT_dds, contrast=c("Tissue", "Thorax", "Ovaries"), alpha = 0.05, lfcThreshold = 1)
summary(thorax_ovaries)
##Order based of padj
ordered_thorax_ovaries <- thorax_ovaries[order(thorax_ovaries$padj),]
##Make data table and write to output
ordered_thorax_ovaries_table <- data.table(data.frame(ordered_thorax_ovaries), keep.rownames = TRUE)
ovaries <- subset(ordered_thorax_ovaries_table, padj < 0.05)

#################################
## overlap for thorax-specific ##
#################################

##venn diagram
vd1 <- venn.diagram(x = list("Head"=head$rn, "Venom"=venom$rn, "Abomen"=abdo$rn, "Ovaries"=ovaries$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd1)

thorax_specific_DEGs <- intersect(intersect(intersect(head$rn, venom$rn), abdo$rn), ovaries$rn)

#####################
## LRT tissue test ##
#####################

mh_dds <- readRDS("output/deseq2/mh_dds.rds")
##select only adult tissue samples
mh_dds_lrt <- mh_dds[,mh_dds$stage=="adult"]
##re-level tissue factor - thorax is tissue of interest so put as first level
mh_dds_lrt$Tissue <- factor(mh_dds_lrt$tissue, levels=c("Thorax", "Head", "Abdomen", "Ovaries", "Venom"))
mh_dds_lrt$Rep <- factor(mh_dds_lrt$rep)
mh_dds_lrt$Flowcell <- factor(mh_dds_lrt$flowcell)
##test
design(mh_dds_lrt) <- ~Flowcell+Rep+Tissue
mh_dds_lrt <- DESeq(mh_dds_lrt, test="LRT", reduced=~Flowcell+Rep)
##results
res_group <- results(mh_dds_lrt, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
lrt_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)

#########################
## overlap WT and LRT  ##
#########################

trinotate <- fread("data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings = ".")
unann_blast_res <- fread("output/deseq2/tissue_itWT_LRT/unann_degs/min_evalues.csv")

LRT_thorax_degs <- subset(lrt_sig_res_group_table, rn %in% thorax_specific_DEGs)
LRT_thorax_degs_annots <- merge(LRT_thorax_degs, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
LRT_thorax_degs_annots_blast <- merge(LRT_thorax_degs_annots, unann_blast_res, by="transcript_id", all.x=TRUE)
fwrite(LRT_thorax_degs_annots_blast, "output/deseq2/tissue_itWT_LRT/thorax/thorax_sp_LRT_annots.csv")

vd_thor <- venn.diagram(x = list("itWT thorax DEGs"=thorax_specific_DEGs, "LRT"=lrt_sig_res_group_table$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_thor)

##not sig in lrt?
not_LRT_thorax_degs_ids <- setdiff(thorax_specific_DEGs$V1, lrt_sig_res_group_table$rn)
thorax_not_lrt <- subset(ordered_res_group_table, rn %in% not_LRT_thorax_degs_ids)

#############
## heatmap ##
#############

##vst transform
mh_vst <- varianceStabilizingTransformation(mh_dds_lrt, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)
##subset for DEGs
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% LRT_thorax_degs$rn)
##turn first row back to row name
mh_vst_degs <- mh_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
##reorder for plot - old samples are off to side for now - need to change this !!!!!!!!!!!!!!!!!!!!!
mh_vst_degs_plot <- mh_vst_degs[,c(10,11,12,4,5,6,16,1,2,3,15,7,8,9,17,18,13,14)]

##get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
##for plot label
sample_to_tissue <- as.data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")])
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
tissue_colours <- list(Tissue = c(Head="#231151FF", Thorax="#5F187FFF", Abdomen="#982D80FF", Ovaries="#D3436EFF", Venom="#F8765CFF"))
##plot
##not clustered by sample
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=FALSE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
##clustered by sample
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
