library(data.table)
library(DESeq2)
library(VennDiagram)

mh_dds <- readRDS("output/deseq2/mh_dds.rds")

##select only adult tissue samples
mh_dds_lrt <- mh_dds[,mh_dds$stage=="adult"]
##re-level tissue factor
mh_dds_lrt$Tissue <- factor(mh_dds_lrt$tissue, levels=c("Head", "Thorax", "Abdomen", "Ovaries", "Venom"))
mh_dds_lrt$Rep <- factor(mh_dds_lrt$rep)
mh_dds_lrt$Flowcell <- factor(mh_dds_lrt$flowcell)

design(mh_dds_lrt) <- ~Flowcell+Rep+Tissue
mh_dds_lrt <- DESeq(mh_dds_lrt, test="LRT", reduced=~Flowcell+Rep)
saveRDS(mh_dds_lrt, "output/deseq2/tissue_LRT/all_mh_tissue_LRT.rds")

mh_dds_lrt <- readRDS("output/deseq2/tissue_LRT/all_mh_tissue_LRT.rds")
res_group <- results(mh_dds_lrt, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_res_group_table, "output/deseq2/tissue_LRT/all_res.csv")
fwrite(ordered_sig_res_group_table, "output/deseq2/tissue_LRT/sig_degs.csv")

plotCounts(mh_dds_lrt, "TRINITY_DN8602_c0_g1", intgroup=("Tissue"))

##compare DEGs found with all samples vs only new samples
all_ordered_sig_res_group_table <- ordered_sig_res_group_table
new_ordered_sig_res_group_table <- fread("output/deseq2/tissue_LRT_only_new/sig_degs.csv")
##venn diagram
vd1 <- venn.diagram(x = list("New"=new_ordered_sig_res_group_table$rn, "All"=all_ordered_sig_res_group_table$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd1)


