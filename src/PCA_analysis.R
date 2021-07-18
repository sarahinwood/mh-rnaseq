library(data.table)
library(DESeq2)
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)

mh_dds <- readRDS("output/deseq2/mh_dds.rds")
mh_dds$tissue <- factor(mh_dds$tissue, levels=c("Head", "Thorax", "Abdomen", "Ovaries", "Venom", "Pupa"))
mh_dds$rep <- factor(mh_dds$rep)

mh_vst <- varianceStabilizingTransformation(mh_dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
pca_plot <- plotPCA(mh_vst, intgroup=c("rep"), returnData=TRUE)
percentVar <- round(100 * attr(pca_plot, "percentVar")) 

##location
##PCA plot (save with dim.s 3.00 x 8.00)
ggplot(pca_plot, aes(x=PC1, y=PC2, color=rep))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis(discrete=TRUE)+
  labs(colour="Rep.s")+
  xlab(paste("PC1:", percentVar[1], "% variance")) + 
  ylab(paste("PC2:", percentVar[2], "% variance")) + 
  coord_fixed()+
  theme_bw()

#########################
## plot cooks distance ##
#########################

mh_dds_LRT <- readRDS("output/deseq2/tissue_LRT/mh_tissue_LRT.rds")
cooks_distance <- data.frame(log10(assays(mh_dds_LRT)[["cooks"]]))
melted_cooks_dist <- cooks_distance %>% gather(colnames(cooks_distance)[1:15], key="sample_name", value="cooks_distance")
sample_table <- fread("data/sample_table.csv")
sample_table$group <- paste(sample_table$tissue)
name_to_group <- sample_table[,c(1, 12)]
cooks_plot <- merge(melted_cooks_dist, name_to_group, all.x=TRUE)
##plot
ggplot(cooks_plot) +
  geom_boxplot(aes(x = sample_name, y = cooks_distance, colour=group), outlier.alpha=0.08) +
  labs(y="Cook's distance", x="", colour="Sample group")+
  scale_colour_viridis(discrete=TRUE)+
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())

boxplot(log10(assays(mh_dds_LRT)[["cooks"]]), range=0, las=2)

##doesn't look like venom3 has consistently higher cook's distance
##some have wider distribution but venom samples narrow

##kraken bacterial results
sample_data <- fread("data/sample_table.csv")
sample_data$tissue <- factor(sample_data$tissue, levels=c("Head", "Thorax", "Abdomen", "Ovaries", "Venom", "Pupa"))
ggplot(sample_data, aes(x=sample_data$tissue, y=sample_data$`kraken_unclassified_%`, colour=tissue))+
  geom_boxplot()+
  labs(x="Tissue", y="Kraken2 unclassified reads (%)")+
  scale_colour_viridis(discrete=TRUE)+
  theme_bw()+
  theme(legend.position = "none")
