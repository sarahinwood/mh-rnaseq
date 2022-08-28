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

library(DESeq2)

###########
# GLOBALS #
###########

mh_dds_file <- snakemake@input[["mh_dds_file"]]

########
# MAIN #
########

mh_dds <- readRDS(mh_dds_file)

##select only adult tissue samples - pupa samples will have all these tissues
mh_dds_tissue <- mh_dds[,mh_dds$stage=="adult"]
##re-level tissue factor
mh_dds_tissue$Tissue <- factor(mh_dds_tissue$tissue, levels=c("Head", "Thorax", "Abdomen", "Ovaries", "Venom"))
mh_dds_tissue$Rep <- factor(mh_dds_tissue$rep)
mh_dds_tissue$Flowcell <- factor(mh_dds_tissue$flowcell)

design(mh_dds_tissue) <- ~Flowcell+Rep+Tissue
mh_dds_tissue <- DESeq(mh_dds_tissue)
saveRDS(mh_dds_tissue, snakemake@output[["mh_itWT_dds"]])

# write log
sessionInfo()
