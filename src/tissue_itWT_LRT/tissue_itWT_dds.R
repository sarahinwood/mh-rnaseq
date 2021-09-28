library(tximport)
library(data.table)
library(DESeq2)
library(EnhancedVolcano)
library(VennDiagram)

mh_dds <- readRDS("output/deseq2/mh_dds.rds")

##select only adult tissue samples - pupa samples will have all these tissues
mh_dds_tissue <- mh_dds[,mh_dds$stage=="adult"]
##re-level tissue factor
mh_dds_tissue$Tissue <- factor(mh_dds_tissue$tissue, levels=c("Head", "Thorax", "Abdomen", "Ovaries", "Venom"))
mh_dds_tissue$Rep <- factor(mh_dds_tissue$rep)
mh_dds_tissue$Flowcell <- factor(mh_dds_tissue$flowcell)

design(mh_dds_tissue) <- ~Flowcell+Rep+Tissue
mh_dds_tissue <- DESeq(mh_dds_tissue)
saveRDS(mh_dds_tissue, "output/deseq2/tissue_itWT_LRT/mh_itWT.rds")

plotCounts(mh_dds_tissue, "TRINITY_DN1993_c0_g1", intgroup=("Tissue"))
