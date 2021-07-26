library(tximport)
library(data.table)
library(DESeq2)
library(EnhancedVolcano)
library(VennDiagram)

mh_dds <- readRDS("output/deseq2/mh_dds.rds")

##select only adult tissue samples - pupa samples will have all these tissues
mh_dds <- mh_dds[,mh_dds$stage=="Adult"]
##re-level tissue factor
mh_dds$Tissue <- factor(mh_dds$tissue, levels=c("Head", "Thorax", "Abdomen", "Ovaries", "Venom"))
mh_dds$Rep <- factor(mh_dds$rep)

design(mh_dds) <- ~Rep+Tissue
mh_dds <- DESeq(mh_dds)
saveRDS(mh_dds, "output/deseq2/itWT/mh_itWT.rds")
