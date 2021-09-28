library(data.table)
library(DESeq2)
library(EnhancedVolcano)
library(VennDiagram)

mh_dds <- readRDS("output/deseq2/mh_dds.rds")
trinotate <- fread("data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings = ".")
viral_annots <- fread("data/mh-transcriptome/output/recip_blast/nr_blastx/viral_annots_plot.csv")

mh_dds$Stage <- factor(mh_dds$stage, levels=c("pupa", "adult"))
mh_dds$Rep <- factor(mh_dds$rep)
mh_dds$Flowcell <- factor(mh_dds$flowcell)

design(mh_dds) <- ~Flowcell+Rep+Stage
mh_dds <- DESeq(mh_dds)
saveRDS(mh_dds, "output/deseq2/stage_WT/mh_stage_WT.rds")

mh_dds <- readRDS("output/deseq2/stage_WT/mh_stage_WT.rds")
res_group <- results(mh_dds, alpha = 0.05, lfcThreshold = 1)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_res_group_table, "output/deseq2/stage_WT/res_group.csv")
fwrite(ordered_sig_res_group_table, "output/deseq2/stage_WT/sig_degs.csv")

EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 1.5, pCutoff = 0.05, colAlpha=0.2, FCcutoff = 1.0,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

##sig DEG annots
sig_annots <- merge(ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id")
fwrite(sig_annots, "output/deseq2/stage_WT/sig_annots.csv")
sig_blastx <- subset(sig_annots, sprot_Top_BLASTX_hit!="")
sig_blastp <- subset(sig_annots, sprot_Top_BLASTP_hit!="")
sig_blast <- full_join(sig_blastx, sig_blastp)
fwrite(sig_blast, "output/deseq2/stage_WT/sig_w_blast_annots.csv")
##any with viral hits in recip blast?
viral_annots$rn <- tstrsplit(viral_annots$transcript_id, "_i", keep=c(1))
viral_degs <- merge(ordered_sig_res_group_table, viral_annots, by="rn")

plotCounts(mh_dds, "TRINITY_DN7695_c0_g1", intgroup=("Stage"))

##tissue specific overlap?
ovary <- fread("output/deseq2/tissue_LRT/ovary/ovary_sp_LRT_annots.csv")
venom <- fread("output/deseq2/tissue_LRT/venom/venom_sp_LRT_annots.csv")
head <- fread("output/deseq2/tissue_LRT/head/head_sp_LRT_annots.csv")
thorax <- fread("output/deseq2/tissue_LRT/thorax/thorax_sp_LRT_annots.csv")
abdo <- fread("output/deseq2/tissue_LRT/abdo/abdo_sp_LRT_annots.csv")

adult <- subset(ordered_sig_res_group_table, log2FoldChange>0)
pupae <- subset(ordered_sig_res_group_table, log2FoldChange<0)

vd_venom <- venn.diagram(x = list("adult"=adult$rn, "pupae"=pupae$rn, "venom"=venom$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_venom)

vd_ovary <- venn.diagram(x = list("adult"=adult$rn, "pupae"=pupae$rn, "ovary"=ovary$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_ovary)

#thorax, abdo, head DEGs not really in pupae and not many in all adult
