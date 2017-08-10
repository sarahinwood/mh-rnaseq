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

##Change what is in these brackets to change which samples you're comparing
#######sting vs head
stingvheadres <- results(dds, contrast = (c("sample", "sting", "head")))
stingvheadresults_table <- data.table(data.frame(stingvheadres), keep.rownames = TRUE)
##me playing around w/sting vs head
stingvhead_resordered <- stingvheadresults_table[order(stingvheadresults_table$padj),]
setkey(stingvhead_resordered, log2FoldChange)
stingvhead_resordered[!is.na(padj)]
fwrite(stingvhead_resordered, "output/deseq2/sting_vs_head_results.csv")
##Generate list of transcripts that are significantly differentially expressed
stingvhead_resSig <- data.table(subset(stingvhead_resordered$rn, stingvhead_resordered$padj < 0.05))
##Save list as .txt
fwrite(stingvhead_resSig, "output/deseq2/stingvhead_sigresults.txt", col.names = FALSE, row.names = FALSE)



#######abdo vs head
abdovheadres <- results(dds, contrast = (c("sample", "abdo", "head")))
abdovhead_results_table <- data.table(data.frame(abdovheadres), keep.rownames = TRUE)
setkey(abdovhead_results_table, log2FoldChange)
abdovhead_results_table[!is.na(padj)]
##
fwrite(abdovhead_results_table, "output/deseq2/abdo_vs_head_results.csv")
##Generate list of transcripts that are significantly differentially expressed
abdovhead_resSig <- as.data.frame(subset(abdovhead_results_table$rn, abdovhead_results_table$padj < 0.05))
##Save list as .txt
fwrite(abdovhead_resSig, "output/deseq2/abdovhead_sigresults.txt", col.names = FALSE, row.names = FALSE)



#######sting vs abdo
res <- results(dds, contrast = (c("sample", "sting", "abdo")))
results_table <- data.table(data.frame(res), keep.rownames = TRUE)
setkey(results_table, log2FoldChange)
results_table[!is.na(padj)]
fwrite(results_table, "output/deseq2/sting_vs_abdo_results.csv")
