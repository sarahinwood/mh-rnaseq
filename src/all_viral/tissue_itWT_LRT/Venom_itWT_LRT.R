library(data.table)
library(DESeq2)
library(pheatmap)
library(VennDiagram)
library(clusterProfiler)

mh_itWT_dds <- readRDS("output/deseq2/viral_LRT/viral_itWT.rds")

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