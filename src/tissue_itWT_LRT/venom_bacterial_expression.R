library(DESeq2)
library(data.table)

venom <- fread("output/03_deseq2/tissue_itWT_LRT/Venom/Venom_sp_LRT_sigp_blast_edited.csv")
venom_bacterial <- subset(venom, venom$`Bacterial?`=="Y")
  
mh_dds_lrt <- readRDS("output/03_deseq2/mh_adult_dds.rds")

#############
## heatmap ##
#############

## vst transform
mh_vst <- varianceStabilizingTransformation(mh_dds_lrt, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)
##subset for DEGs
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% venom_bacterial$rn)
##turn first row back to row name
mh_vst_degs <- mh_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
## reorder tissues for plot - tissue of interest first
mh_vst_degs_plot <- mh_vst_degs[,c(13,14,4,5,6,16,10,11,12,1,2,3,15,7,8,9,17,18)]

## get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mh_dds_lrt)[,c("tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
## for plot label
sample_to_tissue <- as.data.frame(colData(mh_dds_lrt)[,c("tissue", "sample_name")])
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")

tissue_colours <- list(tissue = c(Head='#2D1160FF', Thorax='#721F81FF', Abdomen='#B63679FF', Ovaries='#F1605DFF', Venom='#FEAF77FF'))

## plot
# Modify ordering of the clusters using clustering callback option - pheatmap manual
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
# clustered by sample
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, cluster_cols=TRUE, clustering_callback=callback, show_rownames=FALSE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
