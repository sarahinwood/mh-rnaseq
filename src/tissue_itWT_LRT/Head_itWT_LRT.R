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

##head venom
head_venom <- results(mh_itWT_dds, contrast=c("Tissue", "Head", "Venom"), alpha = 0.05, lfcThreshold = 1)
summary(head_venom)
##Order based of padj
ordered_head_venom <- head_venom[order(head_venom$padj),]
##Make data table and write to output
ordered_head_venom_table <- data.table(data.frame(ordered_head_venom), keep.rownames = TRUE)
venom <- subset(ordered_head_venom_table, padj < 0.05)

##head abdo
head_abdo <- results(mh_itWT_dds, contrast=c("Tissue", "Head", "Abdomen"), alpha = 0.05, lfcThreshold = 1)
summary(head_abdo)
##Order based of padj
ordered_head_abdo <- head_abdo[order(head_abdo$padj),]
##Make data table and write to output
ordered_head_abdo_table <- data.table(data.frame(ordered_head_abdo), keep.rownames = TRUE)
abdo <- subset(ordered_head_abdo_table, padj < 0.05)

##head thorax
head_thorax <- results(mh_itWT_dds, contrast=c("Tissue", "Head", "Thorax"), alpha = 0.05, lfcThreshold = 1)
summary(head_thorax)
##Order based of padj
ordered_head_thorax <- head_thorax[order(head_thorax$padj),]
##Make data table and write to output
ordered_head_thorax_table <- data.table(data.frame(ordered_head_thorax), keep.rownames = TRUE)
thorax <- subset(ordered_head_thorax_table, padj < 0.05)

##head ovaries 
head_ovaries <- results(mh_itWT_dds, contrast=c("Tissue", "Head", "Ovaries"), alpha = 0.05, lfcThreshold = 1)
summary(head_ovaries)
##Order based of padj
ordered_head_ovaries <- head_ovaries[order(head_ovaries$padj),]
##Make data table and write to output
ordered_head_ovaries_table <- data.table(data.frame(ordered_head_ovaries), keep.rownames = TRUE)
ovaries <- subset(ordered_head_ovaries_table, padj < 0.05)

####################################
## itWT overlap for head-specific ##
####################################

##venn diagram
vd1 <- venn.diagram(x = list("Venom"=venom$rn, "Thorax"=thorax$rn, "Abomen"=abdo$rn, "Ovaries"=ovaries$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd1)

head_specific_DEGs <- intersect(intersect(intersect(venom$rn, thorax$rn), abdo$rn), ovaries$rn)

#####################
## LRT tissue test ##
#####################

mh_dds <- readRDS("output/deseq2/mh_dds.rds")
##select only adult tissue samples
mh_dds_lrt <- mh_dds[,mh_dds$stage=="adult"]
##re-level tissue factor - head is tissue of interest so leave as first level
mh_dds_lrt$Tissue <- factor(mh_dds_lrt$tissue, levels=c("Head", "Thorax", "Abdomen", "Ovaries", "Venom"))
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

##annotations
trinotate <- fread("data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings = ".")
unann_blast_res <- fread("output/deseq2/tissue_itWT_LRT/unann_degs/min_evalues.csv")
##subset LRT for head DEGs
LRT_head_degs <- subset(lrt_sig_res_group_table, rn %in% head_specific_DEGs)
##merge with trinotate
LRT_head_degs_annots <- merge(LRT_head_degs, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
##merge with blast
LRT_head_degs_annots_blast <- merge(LRT_head_degs_annots, unann_blast_res, by="transcript_id", all.x=TRUE)
fwrite(LRT_head_degs_annots_blast, "output/deseq2/tissue_itWT_LRT/head/head_sp_LRT_annots.csv")

vd_head <- venn.diagram(x = list("itWT head DEGs"=head_specific_DEGs, "LRT"=lrt_sig_res_group_table$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_head)

##not sig in lrt?
not_LRT_head_degs_ids <- setdiff(head_specific_DEGs, ordered_sig_res_group_table$rn)
head_not_lrt <- subset(ordered_res_group_table, rn %in% not_LRT_head_degs_ids)

#############
## heatmap ##
#############

##vst transform
mh_vst <- varianceStabilizingTransformation(mh_dds_lrt, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)
##subset for DEGs
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% LRT_head_degs$rn)
##turn first row back to row name
mh_vst_degs <- mh_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
##reorder tissues for plot - tissue of interest first
mh_vst_degs_plot <- mh_vst_degs[,c(4,5,6,16,10,11,12,1,2,3,15,7,8,9,17,18,13,14)]

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
