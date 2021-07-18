library(tximport)
library(data.table)
library(DESeq2)
library(EnhancedVolcano)
library(VennDiagram)

mh_dds <- readRDS("output/deseq2/mh_dds.rds")
trinotate <- fread("data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings=".")

##remove outlier venom sample - Mh_venom3
mh_dds_ven <- mh_dds[,-c(18)]

##factors
mh_dds_ven$Tissue <- factor(mh_dds_ven$tissue, levels=c("Head", "Thorax", "Abdomen", "Ovaries", "Venom", "Pupa"))
mh_dds_ven$Rep <- factor(mh_dds_ven$rep)

design(mh_dds_ven) <- ~Rep+Tissue
mh_dds_ven <- DESeq(mh_dds_ven)
#saveRDS(mh_dds, "output/deseq2/venom_itWT/mh_venom_itWT.rds")
#mh_dds <- readRDS("output/deseq2/venom_itWT/mh_venom_itWT.rds")

###############################
## iterative pairwise comp.s ##
###############################

##venom head
venom_head <- results(mh_dds_ven, contrast=c("Tissue", "Venom", "Head"), alpha = 0.05, lfcThreshold = 1)
summary(venom_head)
##Order based of padj
ordered_venom_head <- venom_head[order(venom_head$padj),]
##Make data table and write to output
ordered_venom_head_table <- data.table(data.frame(ordered_venom_head), keep.rownames = TRUE)
head <- subset(ordered_venom_head_table, padj < 0.05)

##venom abdo
venom_abdo <- results(mh_dds_ven, contrast=c("Tissue", "Venom", "Abdomen"), alpha = 0.05, lfcThreshold = 1)
summary(venom_abdo)
##Order based of padj
ordered_venom_abdo <- venom_abdo[order(venom_abdo$padj),]
##Make data table and write to output
ordered_venom_abdo_table <- data.table(data.frame(ordered_venom_abdo), keep.rownames = TRUE)
abdo <- subset(ordered_venom_abdo_table, padj < 0.05)

##venom thorax
venom_thorax <- results(mh_dds_ven, contrast=c("Tissue", "Venom", "Thorax"), alpha = 0.05, lfcThreshold = 1)
summary(venom_thorax)
##Order based of padj
ordered_venom_thorax <- venom_thorax[order(venom_thorax$padj),]
##Make data table and write to output
ordered_venom_thorax_table <- data.table(data.frame(ordered_venom_thorax), keep.rownames = TRUE)
thorax <- subset(ordered_venom_thorax_table, padj < 0.05)

##venom ovaries 
venom_ovaries <- results(mh_dds_ven, contrast=c("Tissue", "Venom", "Ovaries"), alpha = 0.05, lfcThreshold = 1)
summary(venom_ovaries)
##Order based of padj
ordered_venom_ovaries <- venom_ovaries[order(venom_ovaries$padj),]
##Make data table and write to output
ordered_venom_ovaries_table <- data.table(data.frame(ordered_venom_ovaries), keep.rownames = TRUE)
ovaries <- subset(ordered_venom_ovaries_table, padj < 0.05)

##venom pupa 
venom_pupa <- results(mh_dds_ven, contrast=c("Tissue", "Venom", "Pupa"), alpha = 0.05, lfcThreshold = 1)
summary(venom_pupa)
##Order based of padj
ordered_venom_pupa <- venom_pupa[order(venom_pupa$padj),]
##Make data table and write to output
ordered_venom_pupa_table <- data.table(data.frame(ordered_venom_pupa), keep.rownames = TRUE)
pupa <- subset(ordered_venom_pupa_table, padj < 0.05)

################################
## overlap for venom-specific ##
################################

##venn diagram
vd1 <- venn.diagram(x = list("Head"=head$rn, "Thorax"=thorax$rn, "Abomen"=abdo$rn, "Ovaries"=ovaries$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd1)

##1368 venom-specific DEGs without pupa
##when outlier removed - now 733

##venn diagram with pupa
vd2 <- venn.diagram(x = list("Head"=head$rn, "Thorax"=thorax$rn, "Abomen"=abdo$rn, "Ovaries"=ovaries$rn, "Pupa"=pupa$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd2)
##1281 venom-specific DEGs when pupa included - but pupa could contain developing venom glands so don't exclude
##when outlying venom sample is removed from the analysis, down to 704 (lost 577)


venom_specific_DEGs <- intersect(intersect(intersect(head$rn, thorax$rn), abdo$rn), ovaries$rn)
#fwrite(list(venom_specific_DEGs), "output/deseq2/venom_itWT/mh_venom_specific_degs.txt")

venom_sp_annots <- subset(trinotate, `#gene_id` %in% venom_specific_DEGs)

##Still a lot of same annotations for different genes though

