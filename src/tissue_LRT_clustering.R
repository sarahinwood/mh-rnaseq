library(data.table)
library(DESeq2)
library(DEGreport)

##counts and sig genes
tissue_LRT_sig <- fread("output/deseq2/tissue_LRT/sig_degs.csv")
tissue_LRT_rds <- readRDS("output/deseq2/tissue_LRT/mh_tissue_LRT.rds")

##meta data
sample_data <- colData(tissue_LRT_rds)

##need rlog counts for sig genes
rl_dds <- rlog(tissue_LRT_rds, blind=TRUE)
##extract rlog matrix
rld_matrix <- assay(rl_dds)
##filter counts to those for sig DEGs
cluster_rlog <- rld_matrix[tissue_LRT_sig$rn,]

##cluster DEGs
clusters <- degPatterns(cluster_rlog, metadata = sample_data, time = "tissue", col=NULL)

head(clusters$df)
