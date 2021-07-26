library(tximport)
library(data.table)
library(DESeq2)

gene2tx <- fread('data/mh-transcriptome/output/trinity/Trinity.fasta.gene_trans_map', header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files
quant_files <- list.files(path="output/mh_salmon/", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
##assign names to quant files from folder name
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)
##import the salmon quant files (tx2gene links transcript ID to Gene ID
##required for gene-level summarisation for methods that only provide transcript level estimates e.g. salmon)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
##Import table describing samples
sample_data <- fread("data/sample_table.csv", header=TRUE)
setkey(sample_data, sample_name)

##create dds object and link to sample data
dds_all <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
##save dds object
saveRDS(dds_all, 'output/deseq2/mh_dds_all.rds')

##remove venom3 sample
dds_ven <- dds_all[,!(dds_all$sample_name=="Mh_venom3")]
##save dds object
saveRDS(dds_ven, 'output/deseq2/mh_dds.rds')

