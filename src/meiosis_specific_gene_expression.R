library(data.table)
library(tidyverse)
library(pheatmap)
library(viridis)
library(dendsort)

##don't miss any by using length-filtered
trinotate_report <- fread("data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv")
meiosis_genes <- fread("data/meiosis_genes.csv")

##check what genes John looked for and found
meiosis_bx <- subset(trinotate_report, grepl("CDC20|SPO11|HOP2|MND1|DMC1|MSH4|MSH5|REC8|HOP1|RAD21", sprot_Top_BLASTX_hit, ignore.case=TRUE))
meiosis_bp <- subset(trinotate_report, grepl("CDC20|SPO11|HOP2|MND1|DMC1|MSH4|MSH5|REC8|HOP1|RAD21", sprot_Top_BLASTP_hit, ignore.case=TRUE))
trinotate_meiosis <- full_join(meiosis_bx, meiosis_bp)

##remove non-meiosis genes - RMND1 not one of the specifically mentioned genes
trinotate_meiosis_sp <- subset(trinotate_meiosis, !grepl("PHOP1|RMND1", sprot_Top_BLASTX_hit))
trinotate_meiosis_sp$gene <- tstrsplit(trinotate_meiosis_sp$sprot_Top_BLASTX_hit, "_", keep=c(1))
fwrite(trinotate_meiosis_sp, "output/deseq2/tissue_itWT_LRT/ovary/meiosis_genes.csv")

##John found MSH4, MSH5, HOP2 and MND1 conserved at least in microctonus - find all of these
##lrt results for these genes
lrt_degs <- fread("output/deseq2/tissue_itWT_LRT/ovary/sig_ovary_lrt.csv")
meiosis_lrt <- subset(lrt_degs, lrt_degs$rn %in% trinotate_meiosis_sp$`#gene_id`)

##label for DE in LRT
trinotate_meiosis_sp$DE <- ifelse(trinotate_meiosis_sp$`#gene_id` %in% lrt_degs$rn, "*", "")
trinotate_meiosis_sp$Gene <- paste(trinotate_meiosis_sp$gene, trinotate_meiosis_sp$DE, sep="")

gene_id_to_gene <- trinotate_meiosis_sp[,c(2,20)]
setorder(gene_id_to_gene, Gene)
gene_id_to_gene$gene_id_gene <- paste(gene_id_to_gene$Gene, gene_id_to_gene$`#gene_id`, sep="_")
label <- gene_id_to_gene[,c(1,3)]

#############
## heatmap ##
#############

mh_dds_lrt <- readRDS("output/deseq2/tissue_itWT_LRT/ovary/ovary_LRT_dds.rds")
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
##reorder for plot - old samples are off to side for now - need to change this !!!!!!!!!!!!!!!!!!!!!
mh_vst_degs_plot <- mh_vst_degs_label[,c(7,8,9,17,18,4,5,6,16,10,11,12,1,2,3,15,13,14)]

##get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
##for plot label
sample_to_tissue <- as.data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")])
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")

tissue_colours <- list(Tissue = c(Head="#231151FF", Thorax="#5F187FFF", Abdomen="#982D80FF", Ovaries="#D3436EFF", Venom="#F8765CFF"))
##plot
##clustered by sample
col_dend <- dendsort(hclust(dist(t(mh_vst_degs_plot))))
row_dend <- dendsort(hclust(dist(mh_vst_degs_plot)))

pheatmap(mh_vst_degs_plot, cluster_rows=T, cluster_cols=T, show_rownames=TRUE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50), cutree_cols = 2)


#open PDF in inkscape and change gene labels to remove trinity IDs, bold sig DEGs
##dendsort could be used to sort dendrogram? but may be not sorted by longest branch ue to sorting off expression patterns
