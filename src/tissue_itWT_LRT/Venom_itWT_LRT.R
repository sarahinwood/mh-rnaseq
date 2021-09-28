library(data.table)
library(DESeq2)
library(pheatmap)
library(VennDiagram)

mh_itWT_dds <- readRDS("output/deseq2/tissue_itWT_LRT/mh_itWT.rds")

###############################
## iterative pairwise comp.s ##
###############################

##venom head
venom_head <- results(mh_itWT_dds, contrast=c("Tissue", "Venom", "Head"), alpha = 0.05, lfcThreshold = 1)
summary(venom_head)
##Order based of padj
ordered_venom_head <- venom_head[order(venom_head$padj),]
##Make data table and write to output
ordered_venom_head_table <- data.table(data.frame(ordered_venom_head), keep.rownames = TRUE)
head <- subset(ordered_venom_head_table, padj < 0.05)

##venom abdo
venom_abdo <- results(mh_itWT_dds, contrast=c("Tissue", "Venom", "Abdomen"), alpha = 0.05, lfcThreshold = 1)
summary(venom_abdo)
##Order based of padj
ordered_venom_abdo <- venom_abdo[order(venom_abdo$padj),]
##Make data table and write to output
ordered_venom_abdo_table <- data.table(data.frame(ordered_venom_abdo), keep.rownames = TRUE)
abdo <- subset(ordered_venom_abdo_table, padj < 0.05)

##venom thorax
venom_thorax <- results(mh_itWT_dds, contrast=c("Tissue", "Venom", "Thorax"), alpha = 0.05, lfcThreshold = 1)
summary(venom_thorax)
##Order based of padj
ordered_venom_thorax <- venom_thorax[order(venom_thorax$padj),]
##Make data table and write to output
ordered_venom_thorax_table <- data.table(data.frame(ordered_venom_thorax), keep.rownames = TRUE)
thorax <- subset(ordered_venom_thorax_table, padj < 0.05)

##venom ovaries 
venom_ovaries <- results(mh_itWT_dds, contrast=c("Tissue", "Venom", "Ovaries"), alpha = 0.05, lfcThreshold = 1)
summary(venom_ovaries)
##Order based of padj
ordered_venom_ovaries <- venom_ovaries[order(venom_ovaries$padj),]
##Make data table and write to output
ordered_venom_ovaries_table <- data.table(data.frame(ordered_venom_ovaries), keep.rownames = TRUE)
ovaries <- subset(ordered_venom_ovaries_table, padj < 0.05)

################################
## overlap for venom-specific ##
################################

##venn diagram
vd1 <- venn.diagram(x = list("Head"=head$rn, "Thorax"=thorax$rn, "Abomen"=abdo$rn, "Ovaries"=ovaries$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd1)

venom_specific_DEGs <- intersect(intersect(intersect(head$rn, thorax$rn), abdo$rn), ovaries$rn)

#####################
## LRT tissue test ##
#####################

mh_dds <- readRDS("output/deseq2/mh_dds.rds")
##select only adult tissue samples
mh_dds_lrt <- mh_dds[,mh_dds$stage=="adult"]
##re-level tissue factor - thorax is tissue of interest so put as first level
mh_dds_lrt$Tissue <- factor(mh_dds_lrt$tissue, levels=c("Venom", "Head", "Thorax", "Abdomen", "Ovaries"))
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
saveRDS(mh_dds_lrt, "output/deseq2/tissue_itWT_LRT/venom/venom_LRT_dds.rds")


#########################
## overlap WT and LRT  ##
#########################

trinotate <- fread("data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings = ".")
unann_blast_res <- fread("output/deseq2/tissue_itWT_LRT/unann_degs/min_evalues.csv")

LRT_venom_degs <- subset(lrt_sig_res_group_table, rn %in% venom_specific_DEGs)
LRT_venom_degs_annots <- merge(LRT_venom_degs, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
LRT_venom_degs_annots_blast <- merge(LRT_venom_degs_annots, unann_blast_res, by="transcript_id", all.x=TRUE)
fwrite(LRT_venom_degs_annots_blast, "output/deseq2/tissue_itWT_LRT/venom/venom_sp_LRT_annots.csv")

##is signal peptide significantly enriched in venom DEGs?


##annotations seem to be repeated quite a lot?
LRT_venom_degs_annots_blast$blastx <- tstrsplit(LRT_venom_degs_annots_blast$sprot_Top_BLASTX_hit, "^", fixed=TRUE, keep=c(1))
length(unique(LRT_venom_degs_annots_blast$blastx))
blast_annots_count <- count(LRT_venom_degs_annots_blast, blastx)
blast_annots_count <- subset(blast_annots_count, !is.na(blastx))
LRT_venom_degs_annots_blast$blastx_unann <- tstrsplit(LRT_venom_degs_annots_blast$annotation, "<>", fixed=TRUE, keep=c(1))
length(unique(LRT_venom_degs_annots_blast$blastx_unann))
unann_blast_count <- count(LRT_venom_degs_annots_blast, blastx_unann)
unann_blast_count <- subset(unann_blast_count, !is.na(blastx_unann))
mean(unann_blast_count$n)

vd_ven <- venn.diagram(x = list("itWT Venom DEGs"=venom_specific_DEGs, "LRT"=lrt_sig_res_group_table$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_ven)

plotCounts(mh_dds_lrt, "TRINITY_DN25775_c0_g2", intgroup=("Tissue"))

#############
## heatmap ##
#############

mh_dds_lrt <- readRDS("output/deseq2/tissue_itWT_LRT/venom/venom_LRT_dds.rds")
LRT_venom_degs_annots <- fread("output/deseq2/tissue_itWT_LRT/venom/venom_sp_LRT_annots.csv")

##vst transform
mh_vst <- varianceStabilizingTransformation(mh_dds_lrt, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)
##subset for DEGs
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% LRT_venom_degs_annots$rn)
##turn first row back to row name
mh_vst_degs <- mh_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
##reorder for plot
mh_vst_degs_plot <- mh_vst_degs[,c(13,14,4,5,6,16,10,11,12,1,2,3,15,7,8,9,17,18)]

##get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
##for plot label
sample_to_tissue <- as.data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")])
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")

tissue_colours <- list(Tissue = c(Head="#440154FF", Thorax="#3B528BFF", Abdomen="#21908CFF", Ovaries="#5DC863FF", Venom="#FDE725FF"))
##plot
##not clustered by sample
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=FALSE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
##clustered by sample
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
