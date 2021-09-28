library(data.table)
library(DESeq2)

##doesn't find anything meaningful - try stage?

viral_genes_blast <- fread("output/blast/viral_genes/viral_genes_best_hits.csv")
viral_genes_blast$Trinity_gene_ID <- tstrsplit(viral_genes_blast$Trinity_ID, "_i", keep=c(1))
viral_genes <- viral_genes_blast$Trinity_gene_ID

mh_dds <- readRDS("output/deseq2/mh_dds.rds")

##select only adult tissue samples
mh_dds <- mh_dds[,mh_dds$stage=="Adult"]
##keep only viral genes
mh_dds_viral <- mh_dds[viral_genes,]
##re-level tissue factor
mh_dds_viral$Tissue <- factor(mh_dds_viral$tissue, levels=c("Head", "Thorax", "Abdomen", "Ovaries", "Venom"))
mh_dds_viral$Rep <- factor(mh_dds_viral$rep)

design(mh_dds_viral) <- ~Rep+Tissue
##sftype - only works with this setting
mh_dds_viral <- estimateSizeFactors(mh_dds_viral, )
mh_dds_viral <- DESeq(mh_dds_viral, test="LRT", reduced=~Rep, sfType=c("poscounts"))
#saveRDS(mh_dds, "output/deseq2/tissue_LRT/mh_tissue_LRT.rds")

res_group <- results(mh_dds_viral, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/deseq2/tissue_LRT/sig_degs.csv")

plotCounts(mh_dds_viral, "TRINITY_DN25627_c0_g1", intgroup=("Tissue"))

