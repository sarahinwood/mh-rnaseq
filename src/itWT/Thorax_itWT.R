library(tximport)
library(data.table)
library(DESeq2)
library(EnhancedVolcano)
library(VennDiagram)

trinotate <- fread("data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings=".")
mh_dds <- readRDS("output/deseq2/itWT/mh_itWT.rds")

###############################
## iterative pairwise comp.s ##
###############################

##thorax head
thorax_head <- results(mh_dds, contrast=c("Tissue", "Thorax", "Head"), alpha = 0.05, lfcThreshold = 1)
summary(thorax_head)
##Order based of padj
ordered_thorax_head <- thorax_head[order(thorax_head$padj),]
##Make data table and write to output
ordered_thorax_head_table <- data.table(data.frame(ordered_thorax_head), keep.rownames = TRUE)
head <- subset(ordered_thorax_head_table, padj < 0.05)

##thorax abdo
thorax_abdo <- results(mh_dds, contrast=c("Tissue", "Thorax", "Abdomen"), alpha = 0.05, lfcThreshold = 1)
summary(thorax_abdo)
##Order based of padj
ordered_thorax_abdo <- thorax_abdo[order(thorax_abdo$padj),]
##Make data table and write to output
ordered_thorax_abdo_table <- data.table(data.frame(ordered_thorax_abdo), keep.rownames = TRUE)
abdo <- subset(ordered_thorax_abdo_table, padj < 0.05)

##thorax thorax
thorax_venom <- results(mh_dds, contrast=c("Tissue", "Thorax", "Venom"), alpha = 0.05, lfcThreshold = 1)
summary(thorax_venom)
##Order based of padj
ordered_thorax_venom <- thorax_venom[order(thorax_venom$padj),]
##Make data table and write to output
ordered_thorax_venom_table <- data.table(data.frame(ordered_thorax_venom), keep.rownames = TRUE)
venom <- subset(ordered_thorax_venom_table, padj < 0.05)

##thorax ovaries 
thorax_ovaries <- results(mh_dds, contrast=c("Tissue", "Thorax", "Ovaries"), alpha = 0.05, lfcThreshold = 1)
summary(thorax_ovaries)
##Order based of padj
ordered_thorax_ovaries <- thorax_ovaries[order(thorax_ovaries$padj),]
##Make data table and write to output
ordered_thorax_ovaries_table <- data.table(data.frame(ordered_thorax_ovaries), keep.rownames = TRUE)
ovaries <- subset(ordered_thorax_ovaries_table, padj < 0.05)

#################################
## overlap for thorax-specific ##
#################################

##venn diagram
vd1 <- venn.diagram(x = list("Head"=head$rn, "Venom"=venom$rn, "Abomen"=abdo$rn, "Ovaries"=ovaries$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd1)

thorax_specific_DEGs <- intersect(intersect(intersect(head$rn, venom$rn), abdo$rn), ovaries$rn)
fwrite(list(thorax_specific_DEGs), "output/deseq2/itWT/thorax_itWT/mh_thorax_specific_degs.txt")

thorax_sp_annots <- subset(trinotate, `#gene_id` %in% thorax_specific_DEGs)
fwrite(thorax_sp_annots, "output/deseq2/itWT/thorax_itWT/thorax_sp_annots.csv")

