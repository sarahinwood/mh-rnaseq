#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

#############
# LIBRARIES #
#############

library(data.table)
library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(VennDiagram)

###########
# GLOBALS #
###########

mh_dds_file <- snakemake@input[["mh_dds_file"]]
trinotate_file <- snakemake@input[["trinotate_file"]]

########
# MAIN #
########

mh_dds <- readRDS(mh_dds_file)
trinotate <- fread(trinotate_file, na.strings = ".")
#viral_annots <- fread("data/mh-transcriptome/output/recip_blast/nr_blastx/viral_annots_plot.csv")

mh_dds$Stage <- factor(mh_dds$stage, levels=c("pupa", "adult"))
mh_dds$Rep <- factor(mh_dds$rep)
mh_dds$Flowcell <- factor(mh_dds$flowcell)

design(mh_dds) <- ~Flowcell+Rep+Stage
mh_dds <- DESeq(mh_dds)
saveRDS(mh_dds, snakemake@output[["stage_dds"]])

res_group <- results(mh_dds, alpha = 0.05, lfcThreshold = 1)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_res_group_table, snakemake@output[["res_group"]])

pdf(snakemake@output[["volcano"]])
EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 1.5, pCutoff = 0.05, colAlpha=0.2, FCcutoff = 1.0,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))
dev.off()

##sig DEG annots
sig_annots <- merge(ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id")
fwrite(sig_annots, snakemake@output[["sig_annots"]])
sig_blastx <- subset(sig_annots, sprot_Top_BLASTX_hit!="")
sig_blastp <- subset(sig_annots, sprot_Top_BLASTP_hit!="")
sig_blast <- full_join(sig_blastx, sig_blastp)
fwrite(sig_blast, "output/03_deseq2/stage_WT/sig_w_blast_annots.csv")

## any with viral hits in recip blast?
# viral_annots$rn <- tstrsplit(viral_annots$transcript_id, "_i", keep=c(1))
# viral_degs <- merge(ordered_sig_res_group_table, viral_annots, by="rn")

## tissue specific overlap?
ovary <- fread("output/03_deseq2/tissue_itWT_LRT/Ovary/Ovary_sp_LRT_annots.csv")
venom <- fread("output/03_deseq2/tissue_itWT_LRT/Venom/Venom_sp_LRT_annots.csv")
head <- fread("output/03_deseq2/tissue_itWT_LRT/Head/Head_sp_LRT_annots.csv")
thorax <- fread("output/03_deseq2/tissue_itWT_LRT/Thorax/Thorax_sp_LRT_annots.csv")
abdo <- fread("output/03_deseq2/tissue_itWT_LRT/Abdomen/Abdomen_sp_LRT_annots.csv")

adult <- subset(ordered_sig_res_group_table, log2FoldChange>0)
pupae <- subset(ordered_sig_res_group_table, log2FoldChange<0)

pdf(snakemake@output[["venn_venom"]])
vd_venom <- venn.diagram(x = list("adult"=adult$rn, "pupae"=pupae$rn, "venom"=venom$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_venom)
dev.off()

pdf(snakemake@output[["venn_ovary"]])
vd_ovary <- venn.diagram(x = list("adult"=adult$rn, "pupae"=pupae$rn, "ovary"=ovary$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_ovary)
dev.off()

#thorax, abdo, head DEGs not really in pupae and not many in all adult

# write log
sessionInfo()
