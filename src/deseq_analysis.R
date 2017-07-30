library("tximport")
library("data.table")
library("DESeq2")

gene2tx <- fread("data/Trinity.fasta.gene_trans_map", header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

quant_files <- list.files(path="output/salmon", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
sample_data <- data.frame(row.names = names(quant_files), sample=names(quant_files))

dds <- DESeqDataSetFromTximport(txi, colData = sample_data, design = ~sample)
dds <- DESeq(dds)
res <- results(dds)
results_table <- data.table(data.frame(res), keep.rownames = TRUE)
setkey(results_table, log2FoldChange)
results_table[!is.na(padj)]
fwrite(results_table, "output/deseq2/results.csv")
