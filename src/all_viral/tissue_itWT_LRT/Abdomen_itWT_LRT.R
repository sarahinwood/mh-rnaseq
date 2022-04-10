library(data.table)
library(DESeq2)
library(VennDiagram)
library(pheatmap)
library(tidyr)
library(dplyr)
library(viridis)

mh_itWT_dds <- readRDS("output/deseq2/viral_LRT/viral_itWT.rds")

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
