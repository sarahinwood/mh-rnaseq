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
library(VennDiagram)

###########
# GLOBALS #
###########

mh_dds_file <- snakemake@input[["mh_dds_file"]]
trinotate_file <- snakemake@input[["trinotate_file"]]

########
# MAIN #
########

mh_dds <- readRDS(mh_dds_file)

##select only adult tissue samples
mh_dds_lrt <- mh_dds[,mh_dds$stage=="adult"]
##re-level tissue factor
mh_dds_lrt$Tissue <- factor(mh_dds_lrt$tissue, levels=c("Ovaries", "Head", "Thorax", "Abdomen", "Venom"))
mh_dds_lrt$Flowcell <- factor(mh_dds_lrt$flowcell)
mh_dds_lrt$Batch <- factor(mh_dds_lrt$batch)

design(mh_dds_lrt) <- ~Flowcell+Batch+Tissue
mh_dds_lrt <- DESeq(mh_dds_lrt, test="LRT", reduced=~Flowcell+Batch)
saveRDS(mh_dds_lrt, snakemake@output[["dds"]])

res_group <- results(mh_dds_lrt, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_res_group_table, snakemake@output[["all_res"]])
fwrite(ordered_sig_res_group_table, snakemake@output[["sig_res"]])

# write log
sessionInfo()