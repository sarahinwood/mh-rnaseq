library(tximport)
library(data.table)
library(DESeq2)
library(EnhancedVolcano)
library(VennDiagram)

mh_dds <- readRDS("output/deseq2/mh_dds.rds")
trinotate <- fread("data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings=".")

##factors
mh_dds$Tissue <- factor(mh_dds$tissue, levels=c("Head", "Thorax", "Abdomen", "Ovaries", "Venom", "Pupa"))
mh_dds$Rep <- factor(mh_dds$rep)

design(mh_dds) <- ~Rep+Tissue
mh_dds <- DESeq(mh_dds)
saveRDS(mh_dds, "output/deseq2/ovary_itWT/mh_ovaries_itWT.rds")
mh_dds <- readRDS("output/deseq2/ovary_itWT/mh_ovaries_itWT.rds")

###############################
## iterative pairwise comp.s ##
###############################

##venom head
ov_head <- results(mh_dds, contrast=c("Tissue", "Ovaries", "Head"), alpha = 0.05, lfcThreshold = 1)
summary(ov_head)
##Order based of padj
ordered_ov_head <- ov_head[order(ov_head$padj),]
##Make data table and write to output
ordered_ov_head_table <- data.table(data.frame(ordered_ov_head), keep.rownames = TRUE)
head <- subset(ordered_ov_head_table, padj < 0.05)

##venom abdo
ov_abdo <- results(mh_dds, contrast=c("Tissue", "Ovaries", "Abdomen"), alpha = 0.05, lfcThreshold = 1)
summary(ov_abdo)
##Order based of padj
ordered_ov_abdo <- ov_abdo[order(ov_abdo$padj),]
##Make data table and write to output
ordered_ov_abdo_table <- data.table(data.frame(ordered_ov_abdo), keep.rownames = TRUE)
abdo <- subset(ordered_ov_abdo_table, padj < 0.05)

##venom thorax
ov_thorax <- results(mh_dds, contrast=c("Tissue", "Ovaries", "Thorax"), alpha = 0.05, lfcThreshold = 1)
summary(ov_thorax)
##Order based of padj
ordered_ov_thorax <- ov_thorax[order(ov_thorax$padj),]
##Make data table and write to output
ordered_ov_thorax_table <- data.table(data.frame(ordered_ov_thorax), keep.rownames = TRUE)
thorax <- subset(ordered_ov_thorax_table, padj < 0.05)

##venom ovaries 
venom_ov <- results(mh_dds, contrast=c("Tissue", "Ovaries", "Venom"), alpha = 0.05, lfcThreshold = 1)
summary(venom_ov)
##Order based of padj
ordered_venom_ov <- venom_ov[order(venom_ov$padj),]
##Make data table and write to output
ordered_venom_ov_table <- data.table(data.frame(ordered_venom_ov), keep.rownames = TRUE)
venom <- subset(ordered_venom_ov_table, padj < 0.05)

##venom pupa 
ov_pupa <- results(mh_dds, contrast=c("Tissue", "Ovaries", "Pupa"), alpha = 0.05, lfcThreshold = 1)
summary(ov_pupa)
##Order based of padj
ordered_ov_pupa <- ov_pupa[order(ov_pupa$padj),]
##Make data table and write to output
ordered_ov_pupa_table <- data.table(data.frame(ordered_ov_pupa), keep.rownames = TRUE)
pupa <- subset(ordered_ov_pupa_table, padj < 0.05)

################################
## overlap for venom-specific ##
################################

##venn diagram
vd1 <- venn.diagram(x = list("Head"=head$rn, "Thorax"=thorax$rn, "Abomen"=abdo$rn, "Venom"=venom$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd1)

##119 ovary-specific DEGs without pupa

##venn diagram with pupa
vd2 <- venn.diagram(x = list("Head"=head$rn, "Thorax"=thorax$rn, "Abomen"=abdo$rn, "Venom"=venom$rn, "Pupa"=pupa$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd2)
##87 ovary-specific DEGs when pupa included - but pupa could contain developing ovaries so don't exclude

ovary_specific_DEGs <- intersect(intersect(intersect(head$rn, thorax$rn), abdo$rn), venom$rn)
fwrite(list(ovary_specific_DEGs), "output/deseq2/ovary_itWT/mh_ovary_specific_degs.txt")

ovary_sp_annots <- subset(trinotate, `#gene_id` %in% ovary_specific_DEGs)



