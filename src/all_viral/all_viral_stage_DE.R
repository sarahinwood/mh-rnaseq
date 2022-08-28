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

##select only viral genes
mh_dds_wt <- mh_dds[viral_genes$gene_id,]
mh_dds_wt$Stage <- factor(mh_dds_wt$stage)
mh_dds_wt$Rep <- factor(mh_dds_wt$rep)
mh_dds_wt$Flowcell <- factor(mh_dds_wt$flowcell)
##test
design(mh_dds_wt) <- ~Flowcell+Rep+Stage
mh_dds_wt <- DESeq(mh_dds_wt)
##results
res_group <- results(mh_dds_wt, alpha = 0.05, lfcThreshold = 1)
summary(res_group)
saveRDS(mh_dds_wt, "output/03_deseq2/viral_LRT/stage_viral_wt_dds.rds")

##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
wt_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
wt_stage_viral_degs <- merge(wt_sig_res_group_table, viral_genes, by.x="rn", by.y="gene_id", all.x=T)
fwrite(wt_stage_viral_degs, "output/03_deseq2/viral_LRT/stage_viral_degs.csv")

plotCounts(mh_dds_wt, "TRINITY_DN104022_c0_g1", intgroup = ("Stage"))


#############
## heatmap ##
#############

##vst transform
mh_vst <- varianceStabilizingTransformation(mh_dds_wt, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)
##subset for DEGs
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% viral_degs$rn)
##turn first row back to row name
mh_vst_degs <- mh_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")

##reorder for plot
mh_vst_degs_plot <- mh_vst_degs[,c(10:12,1:9,13:21)]

##get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mh_dds_wt)[,c("Stage", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
##for plot label
sample_to_tissue <- as.data.frame(colData(mh_dds_wt)[,c("Stage", "sample_name")])
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
##gene to viral family
gene_to_family <- viral_degs[,c(1,24)]
gene_to_family <- gene_to_family %>% remove_rownames %>% column_to_rownames(var="rn")


tissue_colours <- list(Stage = c(adult="#721F81FF", pupa="#F1605DFF"),
                       family = c(Rhabdoviridae="#0D0887FF", Luteoviridae="#5402A3FF", Nudiviridae="#8B0AA5FF", Baculoviridae="#B93289FF",
                                  Nodaviridae="#DB5C68FF", "Unclassified DNA virus"="#F48849FF", "Unclassified RNA virus"="#FEBC2AFF", Solemoviridae="#F0F921FF"))
##family doesn't match annotation colours - need to fix

##plot
##not clustered by sample
row_dend <- dendsort(hclust(dist(mh_vst_degs_plot)))

pheatmap(mh_vst_degs_plot, cluster_rows=row_dend, cluster_cols=T, show_rownames=F,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_row=gene_to_family,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
