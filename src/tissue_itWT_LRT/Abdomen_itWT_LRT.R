library(data.table)
library(DESeq2)
library(VennDiagram)
library(pheatmap)
library(tidyr)
library(dplyr)
library(viridis)

mh_itWT_dds <- readRDS("output/deseq2/tissue_itWT_LRT/mh_itWT.rds")

###############################
## iterative pairwise comp.s ##
###############################

##abdo head
abdo_head <- results(mh_itWT_dds, contrast=c("Tissue", "Abdomen", "Head"), alpha = 0.05, lfcThreshold = 1)
summary(abdo_head)
##Order based of padj
ordered_abdo_head <- abdo_head[order(abdo_head$padj),]
##Make data table and write to output
ordered_abdo_head_table <- data.table(data.frame(ordered_abdo_head), keep.rownames = TRUE)
head <- subset(ordered_abdo_head_table, padj < 0.05)

##abdo thorax
abdo_thorax <- results(mh_itWT_dds, contrast=c("Tissue", "Abdomen", "Thorax"), alpha = 0.05, lfcThreshold = 1)
summary(abdo_thorax)
##Order based of padj
ordered_abdo_thorax <- abdo_thorax[order(abdo_thorax$padj),]
##Make data table and write to output
ordered_abdo_thorax_table <- data.table(data.frame(ordered_abdo_thorax), keep.rownames = TRUE)
thorax <- subset(ordered_abdo_thorax_table, padj < 0.05)

##abdo thorax
abdo_venom <- results(mh_itWT_dds, contrast=c("Tissue", "Abdomen", "Venom"), alpha = 0.05, lfcThreshold = 1)
summary(abdo_venom)
##Order based of padj
ordered_abdo_venom <- abdo_venom[order(abdo_venom$padj),]
##Make data table and write to output
ordered_abdo_venom_table <- data.table(data.frame(ordered_abdo_venom), keep.rownames = TRUE)
venom <- subset(ordered_abdo_venom_table, padj < 0.05)

##abdo ovaries 
abdo_ovaries <- results(mh_itWT_dds, contrast=c("Tissue", "Abdomen", "Ovaries"), alpha = 0.05, lfcThreshold = 1)
summary(abdo_ovaries)
##Order based of padj
ordered_abdo_ovaries <- abdo_ovaries[order(abdo_ovaries$padj),]
##Make data table and write to output
ordered_abdo_ovaries_table <- data.table(data.frame(ordered_abdo_ovaries), keep.rownames = TRUE)
ovaries <- subset(ordered_abdo_ovaries_table, padj < 0.05)

##################################
## overlap for abdomen-specific ##
##################################

##venn diagram
## abdomen-specific DEGs
vd1 <- venn.diagram(x = list("Head"=head$rn, "Venom"=venom$rn, "Thorax"=thorax$rn, "Ovaries"=ovaries$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd1)

abdo_specific_DEGs <- intersect(intersect(intersect(head$rn, venom$rn), thorax$rn), ovaries$rn)

#####################
## LRT tissue test ##
#####################

mh_dds <- readRDS("output/deseq2/mh_dds.rds")
##select only adult tissue samples
mh_dds_lrt <- mh_dds[,mh_dds$stage=="adult"]
##re-level tissue factor - thorax is tissue of interest so put as first level
mh_dds_lrt$Tissue <- factor(mh_dds_lrt$tissue, levels=c("Abdomen","Head", "Thorax", "Ovaries", "Venom"))
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
LRT_abdo_degs <- subset(lrt_sig_res_group_table, rn %in% abdo_specific_DEGs)
LRT_abdo_degs_annots <- merge(LRT_abdo_degs, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
LRT_abdo_degs_annots_blast <- merge(LRT_abdo_degs_annots, unann_blast_res, by="transcript_id", all.x=TRUE)
fwrite(LRT_abdo_degs_annots_blast, "output/deseq2/tissue_itWT_LRT/abdo/abdo_sp_LRT_annots.csv")

vd_abdo <- venn.diagram(x = list("itWT abdo DEGs"=abdo_specific_DEGs, "LRT"=lrt_sig_res_group_table$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_abdo)

#############
## heatmap ##
#############

##vst transform
mh_vst <- varianceStabilizingTransformation(mh_dds_lrt, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)
##subset for DEGs
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% abdo_specific_DEGs)
##turn first row back to row name
mh_vst_degs <- mh_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
##reorder for plot - old samples are off to side for now - need to change this !!!!!!!!!!!!!!!!!!!!!
mh_vst_degs_plot <- mh_vst_degs[,c(1,2,3,4,16,5,6,7,17,11,12,13,8,9,10,18,19,14,15)]
mh_vst_degs_plot <- mh_vst_degs_plot %>% column_to_rownames(var="rn")

##get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
##for plot label
sample_to_tissue <- as.data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")])
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")

##why is thos not working?
tissue_colours <- list(Tissue = c(Head="#440154FF", Thorax="#3B528BFF", Abdomen="#21908CFF", Ovaries="#5DC863FF", Venom="#FDE725FF"))
##plot
##not clustered by samples
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=FALSE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
##clustered by samples
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
##clustered by samples with gene names
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))

##heatmap.2 instead - similar result
library(gplots)
heatmap.2(as.matrix(mh_vst_degs_plot))


