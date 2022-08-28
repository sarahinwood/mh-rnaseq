library(data.table)
library(tidyverse)
library(pheatmap)
library(viridis)

mh_dds <- readRDS("output/03_deseq2/mh_dds.rds")
mv_dw_m_virus <- fread("data/mh-all-transcriptome/output/recip_blast/viral_nr_blastx/iflavirus_blastn_best.csv")
mv_dw_m_virus$gene_id <- tstrsplit(mv_dw_m_virus$transcript_id, "_i", keep=c(1))

#############
## heatmap ##
#############

##vst transform
mh_vst <- varianceStabilizingTransformation(mh_dds, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)
##subset for DEGs
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% mv_dw_m_virus$gene_id)
##turn first row back to row name
mh_vst_degs <- mh_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")

##reorder for plot
mh_vst_degs_plot <- mh_vst_degs[,c(4,5,6,19,13,14,15,1,2,3,18,7,8,9,20,21,16,17,10,11,12)]

##get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mh_dds)[,c("tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
##for plot label
sample_to_tissue <- as.data.frame(colData(mh_dds)[,c("tissue", "sample_name")])
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")

tissue_colours <- list(tissue = c(Head="#231151FF", Thorax="#5F187FFF", Abdomen="#982D80FF", Ovaries="#D3436EFF", Venom="#F8765CFF", Pupa="#FCFDBFFF"))
##plot
pheatmap(mh_vst_degs_plot, cluster_rows=T, cluster_cols=F, show_rownames=T,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
