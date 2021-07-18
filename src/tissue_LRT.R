library(tximport)
library(data.table)
library(DESeq2)
library(VennDiagram)

mh_dds <- readRDS("output/deseq2/mh_dds.rds")
trinotate <- fread("data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings = ".")

##select only adult tissue samples
mh_dds <- mh_dds[,mh_dds$stage=="Adult"]
##re-level tissue factor
mh_dds$tissue <- factor(mh_dds$tissue, levels=c("Head", "Thorax", "Abdomen", "Ovaries", "Venom"))
mh_dds$Rep <- factor(mh_dds$rep)

design(mh_dds) <- ~Rep+Tissue
mh_dds <- DESeq(mh_dds, test="LRT", reduced=~Rep)
saveRDS(mh_dds, "output/deseq2/tissue_LRT/mh_tissue_LRT.rds")

res_group <- results(mh_dds, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/deseq2/tissue_LRT/sig_degs.csv")

##gives error about it being a large number of genes to work with
plotCounts(mh_dds, "TRINITY_DN5452_c0_g1", intgroup=("tissue"))

################################
## tissue-specific gene lists ##
################################

##ovary
ovary_specific_DEGs <- fread("output/deseq2/ovary_itWT/mh_ovary_specific_degs.txt", header=FALSE)
LRT_ovary_degs <- subset(ordered_sig_res_group_table, rn %in% ovary_specific_DEGs$V1)
LRT_ovary_degs_annots <- merge(LRT_ovary_degs, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
##87 DEGs
vd_ov <- venn.diagram(x = list("itWT Ovary DEGs"=ovary_specific_DEGs$V1, "LRT"=ordered_sig_res_group_table$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_ov)
##all itWT DEGs found by LRT

##venom
venom_specific_DEGs <- fread("output/deseq2/venom_itWT/mh_venom_specific_degs.txt", header=FALSE)
LRT_venom_degs <- subset(ordered_sig_res_group_table, rn %in% venom_specific_DEGs$V1)
LRT_venom_degs_annots <- merge(LRT_venom_degs, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
##1181 DEGs
vd_ven <- venn.diagram(x = list("itWT Venom DEGs"=venom_specific_DEGs$V1, "LRT"=ordered_sig_res_group_table$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_ven)
##187 DEGs found from itWT not found in LRT - why??
itWT_ven_sp <- setdiff(venom_specific_DEGs$V1, ordered_sig_res_group_table$rn)
itWT_ven_res <- subset(ordered_res_group_table, rn %in% itWT_ven_sp)
