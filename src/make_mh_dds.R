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

library(tximport)
library(data.table)
library(DESeq2)

###########
# GLOBALS #
###########

gene2tx_file <- snakemake@input[["gene2tx_file"]]
sample_data_file <- snakemake@input[["sample_data_file"]]

########
# MAIN #
########

gene2tx <- fread(gene2tx_file, header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files
quant_files <- list.files(path="output/01_mh_salmon/", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
##assign names to quant files from folder name
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)
##import the salmon quant files (tx2gene links transcript ID to Gene ID
##required for gene-level summarisation for methods that only provide transcript level estimates e.g. salmon)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE, abundanceCol="TPM")
txi_abundance <- data.table(txi$abundance, keep.rownames=TRUE)
fwrite(txi_abundance, snakemake@output[["salmon_tpm"]])

##Import table describing samples
sample_data <- fread(sample_data_file, header=TRUE)
setkey(sample_data, sample_name)

##create dds object and link to sample data
mh_dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
##save dds object
saveRDS(mh_dds, snakemake@output[["mh_dds"]])

mh_adult_dds <- mh_dds_lrt <- mh_dds[,mh_dds$stage=="adult"]
saveRDS(mh_adult_dds, snakemake@output[["mh_adult_dds"]])

# write log
sessionInfo()

