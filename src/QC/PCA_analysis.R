#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

#############
# LIBRARIES #
#############

library(data.table)
library(DESeq2)
library(ggplot2)
library(viridis)
library(tidyverse)

###########
# GLOBALS #
###########

mh_dds_file <- snakemake@input[["mh_dds_file"]]
sample_data_file <- snakemake@input[["sample_data_file"]]

########
# MAIN #
########

mh_dds <- readRDS(mh_dds_file)
sample_data <- fread(sample_data_file)
sample_to_seqrun <- sample_data[,c(1,5)]

mh_dds$tissue <- factor(mh_dds$tissue, levels=c("Head", "Thorax", "Abdomen", "Ovaries", "Venom", "Pupa"))
mh_dds$rep <- factor(mh_dds$rep)
mh_dds$flowcell <- factor(mh_dds$flowcell)

mh_vst <- varianceStabilizingTransformation(mh_dds, blind=TRUE)

##plot PCA with first 2 dimensions to investigate sample clustering - tissue
pca_plot_tissue <- plotPCA(mh_vst, intgroup=c("tissue"), returnData=TRUE)
percentVar <- round(100 * attr(pca_plot_tissue, "percentVar")) 
##PCA plot (save with dim.s 3.00 x 8.00)
## tissue ##
pdf(snakemake@output[["PCA_tissue"]])
ggplot(pca_plot_tissue, aes(x=PC1, y=PC2, color=tissue))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis(discrete=TRUE)+
  labs(colour="Tissue")+
  xlab(paste("PC1:", percentVar[1], "% variance")) + 
  ylab(paste("PC2:", percentVar[2], "% variance")) + 
  coord_fixed()+
  theme_bw()
dev.off()

## replicate ##
pca_plot_rep <- plotPCA(mh_vst, intgroup=c("rep"), returnData=TRUE)
percentVar <- round(100 * attr(pca_plot_rep, "percentVar")) 
pdf(snakemake@output[["PCA_replicate"]])
ggplot(pca_plot_rep, aes(x=PC1, y=PC2, color=rep))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis(discrete=TRUE)+
  labs(colour="Replicate")+
  xlab(paste("PC1:", percentVar[1], "% variance")) + 
  ylab(paste("PC2:", percentVar[2], "% variance")) + 
  coord_fixed()+
  theme_bw()
dev.off()

## flowcell ##
pca_plot_run <- plotPCA(mh_vst, intgroup=c("flowcell"), returnData=TRUE)
percentVar <- round(100 * attr(pca_plot_run, "percentVar"))
pdf(snakemake@output[["PCA_seqrun"]])
ggplot(pca_plot_run, aes(x=PC1, y=PC2, color=flowcell))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis(discrete=TRUE)+
  labs(colour="Sequencing run")+
  xlab(paste("PC1:", percentVar[1], "% variance")) + 
  ylab(paste("PC2:", percentVar[2], "% variance")) + 
  coord_fixed()+
  theme_bw()
dev.off()

# write log
sessionInfo()
