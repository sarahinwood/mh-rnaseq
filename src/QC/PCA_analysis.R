library(data.table)
library(DESeq2)
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)

##can use on all samples or with venom3 removed
mh_dds <- readRDS("output/deseq2/mh_dds_all.rds")
mh_dds <- readRDS("output/deseq2/mh_dds.rds")

mh_dds$tissue <- factor(mh_dds$tissue, levels=c("Head", "Thorax", "Abdomen", "Ovaries", "Venom", "Pupa"))
mh_dds$rep <- factor(mh_dds$rep)

mh_vst <- varianceStabilizingTransformation(mh_dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
pca_plot <- plotPCA(mh_vst, intgroup=c("tissue"), returnData=TRUE)
percentVar <- round(100 * attr(pca_plot, "percentVar")) 

##PCA plot (save with dim.s 3.00 x 8.00)
ggplot(pca_plot, aes(x=PC1, y=PC2, color=tissue))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis(discrete=TRUE)+
  labs(colour="Tissue")+
  xlab(paste("PC1:", percentVar[1], "% variance")) + 
  ylab(paste("PC2:", percentVar[2], "% variance")) + 
  coord_fixed()+
  theme_bw()
