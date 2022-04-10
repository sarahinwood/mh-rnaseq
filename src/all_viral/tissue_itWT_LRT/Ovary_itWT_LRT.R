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

plotCounts(mh_itWT_dds, "TRINITY_DN44471_c0_g2", intgroup = ("Tissue"))
