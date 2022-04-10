library(data.table)
library(DESeq2)
library(VennDiagram)
library(pheatmap)
library(viridis)

mh_itWT_dds <- readRDS("output/deseq2/tissue_itWT_LRT/mh_itWT.rds")

###############################
## iterative pairwise comp.s ##
###############################

##ovaries head
ov_head <- results(mh_itWT_dds, contrast=c("Tissue", "Ovaries", "Head"), alpha = 0.05, lfcThreshold = 1)
summary(ov_head)
##Order based of padj
ordered_ov_head <- ov_head[order(ov_head$padj),]
##Make data table and write to output
ordered_ov_head_table <- data.table(data.frame(ordered_ov_head), keep.rownames = TRUE)
head <- subset(ordered_ov_head_table, padj < 0.05)

##ovaries abdo
ov_abdo <- results(mh_itWT_dds, contrast=c("Tissue", "Ovaries", "Abdomen"), alpha = 0.05, lfcThreshold = 1)
summary(ov_abdo)
##Order based of padj
ordered_ov_abdo <- ov_abdo[order(ov_abdo$padj),]
##Make data table and write to output
ordered_ov_abdo_table <- data.table(data.frame(ordered_ov_abdo), keep.rownames = TRUE)
abdo <- subset(ordered_ov_abdo_table, padj < 0.05)

##ovaries thorax
ov_thorax <- results(mh_itWT_dds, contrast=c("Tissue", "Ovaries", "Thorax"), alpha = 0.05, lfcThreshold = 1)
summary(ov_thorax)
##Order based of padj
ordered_ov_thorax <- ov_thorax[order(ov_thorax$padj),]
##Make data table and write to output
ordered_ov_thorax_table <- data.table(data.frame(ordered_ov_thorax), keep.rownames = TRUE)
thorax <- subset(ordered_ov_thorax_table, padj < 0.05)

##ovaries venom 
ov_venom <- results(mh_itWT_dds, contrast=c("Tissue", "Ovaries", "Venom"), alpha = 0.05, lfcThreshold = 1)
summary(ov_venom)
##Order based of padj
ordered_venom_ov <- ov_venom[order(ov_venom$padj),]
##Make data table and write to output
ordered_venom_ov_table <- data.table(data.frame(ordered_venom_ov), keep.rownames = TRUE)
venom <- subset(ordered_venom_ov_table, padj < 0.05)

################################
## overlap for ovary-specific ##
################################

##venn diagram
vd1 <- venn.diagram(x = list("Head"=head$rn, "Thorax"=thorax$rn, "Abomen"=abdo$rn, "Venom"=venom$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd1)

ovary_specific_DEGs <- intersect(intersect(intersect(head$rn, thorax$rn), abdo$rn), venom$rn)

#####################
## LRT tissue test ##
#####################

mh_dds <- readRDS("output/deseq2/mh_dds.rds")
##select only adult tissue samples
mh_dds_lrt <- mh_dds[,mh_dds$stage=="adult"]
##re-level tissue factor - thorax is tissue of interest so put as first level
mh_dds_lrt$Tissue <- factor(mh_dds_lrt$tissue, levels=c("Ovaries", "Head", "Thorax", "Abdomen", "Venom"))
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
fwrite(lrt_sig_res_group_table, "output/deseq2/tissue_itWT_LRT/ovary/sig_ovary_lrt.csv")
saveRDS(mh_dds_lrt, "output/deseq2/tissue_itWT_LRT/ovary/ovary_LRT_dds.rds")

#########################
## overlap WT and LRT  ##
#########################

trinotate <- fread("data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings = ".")
unann_blast_res <- fread("output/deseq2/tissue_itWT_LRT/unann_degs/min_evalues.csv")

LRT_ovary_degs <- subset(lrt_sig_res_group_table, rn %in% ovary_specific_DEGs)
LRT_ovary_degs_annots <- merge(LRT_ovary_degs, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
LRT_ovary_degs_annots_blast <- merge(LRT_ovary_degs_annots, unann_blast_res, by="transcript_id", all.x=TRUE)
fwrite(LRT_ovary_degs_annots_blast, "output/deseq2/tissue_itWT_LRT/ovary/ovary_sp_LRT_annots.csv")

vd_ov <- venn.diagram(x = list("itWT Ovary DEGs"=ovary_specific_DEGs, "LRT"=lrt_sig_res_group_table$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_ov)

#############
## heatmap ##
#############

##vst transform
mh_vst <- varianceStabilizingTransformation(mh_dds_lrt, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)
##subset for DEGs
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% ovary_specific_DEGs)
##turn first row back to row name
mh_vst_degs <- mh_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
##reorder for plot - old samples are off to side for now - need to change this !!!!!!!!!!!!!!!!!!!!!
mh_vst_degs_plot <- mh_vst_degs[,c(7,8,9,17,18,4,5,6,16,10,11,12,1,2,3,15,13,14)]

##get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
##for plot label
sample_to_tissue <- as.data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")])
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")

tissue_colours <- list(Tissue = c(Head="#440154FF", Thorax="#3B528BFF", Abdomen="#21908CFF", Ovaries="#5DC863FF", Venom="#FDE725FF"), cluster=c("1"="#000004FF", "2"="#231151FF", "3"="#5F187FFF", "4"="#982D80FF", "5"="#D3436EFF", "6"="#F8765CFF", "7"="#FEBA80FF", "8"="#FCFDBFFF"))
##plot
##not clustered by sample
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=FALSE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
##clustered by sample
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
#cutting off abdomen label from graph - open PDF in inkscape and resize page

##write table of geneID to clusters -make call on what k= for clusters from plot
pheatmap_dt <- pheatmap(mh_vst_degs_plot)
pheatmap_clusters <- data.table(cbind(mh_vst_degs_plot, cluster=cutree(pheatmap_dt$tree_row, k=8)), keep.rownames = TRUE)
id_to_cluster <- pheatmap_clusters[,c(1,20)]

id_to_cluster <- id_to_cluster %>% column_to_rownames(var="rn")

fwrite(id_to_cluster, "output/deseq2/tissue_itWT_LRT/ovary/ovary_DEG_clusters.csv")
fwrite(pheatmap_clusters,)

vst_degs_dt <- data.table(mh_vst_degs_plot, keep.rownames = TRUE)
vst_clusters <- merge(vst_degs_dt, id_to_cluster, by="rn")
vst_plot <- vst_clusters %>% column_to_rownames(var="cluster")
##clustered by sample
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50), cutree_rows=8, annotation_row=id_to_cluster)


