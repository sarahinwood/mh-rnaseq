library(DESeq2)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(viridis)

mh_dds_lrt <- readRDS("output/deseq2/tissue_LRT/all_mh_tissue_LRT.rds")
##read in degs for specific tissue
tissue_degs <- fread("output/deseq2/tissue_LRT/head/head_sp_LRT_annots.csv")

##vst transform
mh_vst <- varianceStabilizingTransformation(mh_dds_lrt, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)
##subset for DEGs
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% tissue_degs$rn)
##turn first row back to row name
mh_vst_degs <- mh_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
##reorder for plot - old samples are off to side for now
mh_vst_degs_plot <- mh_vst_degs[,c(4,5,6,16,10,11,12,1,2,3,15,7,8,9,17,18,13,14)]

##get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
##for plot label
sample_to_tissue <- as.data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")])
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")

##plot
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=sample_to_tissue, color=viridis(50))


