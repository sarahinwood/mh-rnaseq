library(tximport)
library(data.table)
library(DESeq2)
library(EnhancedVolcano)

mh_dds <- readRDS("output/deseq2/mh_dds.rds")
trinotate <- fread("data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings = ".")

mh_dds$Stage <- factor(mh_dds$stage, levels=c("Pupa", "Adult"))
mh_dds$Rep <- factor(mh_dds$rep)

design(mh_dds) <- ~Rep+Stage
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

plotCounts(mh_dds, "TRINITY_DN54112_c4_g1", intgroup=("Stage"))

EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 1.5, pCutoff = 0.05, colAlpha=0.2, FCcutoff = 1.0,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

##sig DEG annots
sig_annots <- merge(ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id")
fwrite(sig_annots, "output/deseq2/stage_WT/sig_annots.csv")

