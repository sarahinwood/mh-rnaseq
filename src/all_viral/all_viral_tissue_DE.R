library(data.table)
library(tidyverse)
library(pheatmap)
library(viridis)
library(dendsort)

viral_genes <- fread("data/mh-all-transcriptome/output/recip_blast/viral_nr_blastx/best_viral_hits_plot.csv")
viral_genes$gene_id <- tstrsplit(viral_genes$transcript_id, "_i", keep=c(1))

mh_dds <- readRDS("output/03_deseq2/mh_dds.rds")

########
## DE ##
########

##only adult for tissue analysis
mh_dds_lrt <- mh_dds[,mh_dds$stage=="adult"]
##select only viral genes
mh_dds_lrt <- mh_dds_lrt[viral_genes$gene_id,]
mh_dds_lrt$Tissue <- factor(mh_dds_lrt$tissue)
mh_dds_lrt$Rep <- factor(mh_dds_lrt$rep)
mh_dds_lrt$Flowcell <- factor(mh_dds_lrt$flowcell)
##test
design(mh_dds_lrt) <- ~Flowcell+Rep+Tissue
mh_dds_lrt <- DESeq(mh_dds_lrt, test="LRT", reduced=~Flowcell+Rep)
##results
res_group <- results(mh_dds_lrt, alpha = 0.05)
summary(res_group)
saveRDS(mh_dds_lrt, "output/03_deseq2/viral_LRT/viral_LRT_dds.rds")

##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
lrt_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
viral_degs <- merge(lrt_sig_res_group_table, viral_genes, by.x="rn", by.y="gene_id", all.x=T)
fwrite(viral_degs, "output/03_deseq2/viral_LRT/tissue_lrt_viral_degs.csv")

#############
## heatmap ##
#############

##vst transform
mh_vst <- varianceStabilizingTransformation(mh_dds_lrt, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)
##subset for DEGs
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% viral_degs$rn)
##turn first row back to row name
mh_vst_degs <- mh_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")

##reorder for plot
mh_vst_degs_plot <- mh_vst_degs[,c(4,5,6,16,10,11,12,1,2,3,15,7,8,9,17,18,13,14)]

##get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
##for plot label
sample_to_tissue <- as.data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")])
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
##gene to viral family
gene_to_family <- viral_degs[,c(1,24)]
gene_to_family <- gene_to_family %>% remove_rownames %>% column_to_rownames(var="rn")


tissue_colours <- list(Tissue = c(Head="#231151FF", Thorax="#5F187FFF", Abdomen="#982D80FF", Ovaries="#D3436EFF", Venom="#F8765CFF"),
                       family = c(Baculoviridae="#0D0887FF", Iflaviridae="#5D01A6FF", Nudiviridae="#9C179EFF", Polydnaviridae="#CC4678FF",
                                  Rhabdoviridae="#ED7953FF", "Unclassified DNA virus"="#FDB32FFF", "Unclassified RNA virus"="#F0F921FF"))
##plot
##not clustered by sample
row_dend <- dendsort(hclust(dist(mh_vst_degs_plot)))

pheatmap(mh_vst_degs_plot, cluster_rows=row_dend, cluster_cols=F, show_rownames=F,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_row=gene_to_family,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
