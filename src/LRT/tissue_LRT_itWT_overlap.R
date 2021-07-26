library(tximport)
library(data.table)
library(DESeq2)
library(VennDiagram)

trinotate <- fread("data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings = ".")
##LRT results
ordered_sig_res_group_table <- fread("output/deseq2/tissue_LRT/sig_degs.csv")
##reduce to 0.01 and see if overlap is maintained
ordered_sig_res_group_table_01 <- subset(ordered_sig_res_group_table, padj<0.01)

##compare genes found in tissue itWT with LRT DEGs

###########
## ovary ##
###########
ovary_specific_DEGs <- fread("output/deseq2/itWT/ovary_itWT/mh_ovary_specific_degs.txt", header=FALSE)
vd_ov <- venn.diagram(x = list("itWT Ovary DEGs"=ovary_specific_DEGs$V1, "LRT"=ordered_sig_res_group_table_01$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_ov)
##all itWT DEGs found by LRT
LRT_ovary_degs <- subset(ordered_sig_res_group_table, rn %in% ovary_specific_DEGs$V1)
LRT_ovary_degs_annots <- merge(LRT_ovary_degs, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(LRT_ovary_degs_annots, "output/deseq2/tissue_LRT/ovary/ovary_sp_LRT_annots.csv")

###########
## venom ##
###########
venom_specific_DEGs <- fread("output/deseq2/itWT/venom_itWT/mh_venom_specific_degs.txt", header=FALSE)
vd_ven <- venn.diagram(x = list("itWT Venom DEGs"=venom_specific_DEGs$V1, "LRT"=ordered_sig_res_group_table_01$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_ven)
##all but 1 DEGs found by LRT
LRT_venom_degs <- subset(ordered_sig_res_group_table, rn %in% venom_specific_DEGs$V1)
LRT_venom_degs_annots <- merge(LRT_venom_degs, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(LRT_venom_degs_annots, "output/deseq2/tissue_LRT/venom/venom_sp_LRT_annots.csv")

##########
## head ##
##########
head_specific_DEGs <- fread("output/deseq2/itWT/head_itWT/mh_head_specific_degs.txt", header=FALSE)
vd_head <- venn.diagram(x = list("itWT head DEGs"=head_specific_DEGs$V1, "LRT"=ordered_sig_res_group_table_01$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_head)
LRT_head_degs <- subset(ordered_sig_res_group_table, rn %in% head_specific_DEGs$V1)
LRT_head_degs_annots <- merge(LRT_head_degs, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(LRT_head_degs_annots, "output/deseq2/tissue_LRT/head/head_sp_LRT_annots.csv")

############
## thorax ##
############
thorax_specific_DEGs <- fread("output/deseq2/itWT/thorax_itWT/mh_thorax_specific_degs.txt", header=FALSE)
vd_thor <- venn.diagram(x = list("itWT thorax DEGs"=thorax_specific_DEGs$V1, "LRT"=ordered_sig_res_group_table_01$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_thor)
LRT_thorax_degs <- subset(ordered_sig_res_group_table, rn %in% thorax_specific_DEGs$V1)
LRT_thorax_degs_annots <- merge(LRT_thorax_degs, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(LRT_thorax_degs_annots, "output/deseq2/tissue_LRT/thorax/thorax_sp_LRT_annots.csv")

##########
## abdo ##
##########
abdo_specific_DEGs <- fread("output/deseq2/itWT/abdo_itWT/mh_abdo_specific_degs.txt", header=FALSE)
vd_abdo <- venn.diagram(x = list("itWT abdo DEGs"=abdo_specific_DEGs$V1, "LRT"=ordered_sig_res_group_table_01$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_abdo)
LRT_abdo_degs <- subset(ordered_sig_res_group_table, rn %in% abdo_specific_DEGs$V1)
LRT_abdo_degs_annots <- merge(LRT_abdo_degs, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(LRT_abdo_degs_annots, "output/deseq2/tissue_LRT/abdo/abdo_sp_LRT_annots.csv")

##overlap to find what number of DEGs from LRT are NOT found in tissue iTWT?
