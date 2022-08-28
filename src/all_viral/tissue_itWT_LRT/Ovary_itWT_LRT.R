library(data.table)
library(DESeq2)
library(VennDiagram)
library(pheatmap)
library(viridis)

viral_genes <- fread("data/mh-all-transcriptome/output/recip_blast/viral_nr_blastx/best_viral_hits_plot.csv")
viral_genes$gene_id <- tstrsplit(viral_genes$transcript_id, "_i", keep=c(1))
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

mh_dds_lrt <- readRDS("output/03_deseq2/viral_LRT/viral_LRT_dds.rds")
## re-level tissue factor - ovaries are tissue of interest so put as first level for lfc calculation
mh_dds_lrt$Tissue <- factor(mh_dds_lrt$tissue, levels=c("Ovaries", "Head", "Thorax", "Abdomen", "Venom"))
mh_dds_lrt$Rep <- factor(mh_dds_lrt$rep)
mh_dds_lrt$Flowcell <- factor(mh_dds_lrt$flowcell)
## test
design(mh_dds_lrt) <- ~Flowcell+Rep+Tissue
mh_dds_lrt <- DESeq(mh_dds_lrt, test="LRT", reduced=~Flowcell+Rep)
## results
res_group <- results(mh_dds_lrt, alpha = 0.05)
summary(res_group)
## Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
## Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
lrt_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
#### overlap WT and LRT
LRT_WT_ovary_degs <- subset(lrt_sig_res_group_table, rn %in% ovary_specific_DEGs)
LRT_WT_ovary_degs_annots <- merge(LRT_WT_ovary_degs, viral_genes, by.x="rn", by.y="gene_id")

plotCounts(mh_itWT_dds, "TRINITY_DN44471_c0_g2", intgroup = ("Tissue"))


