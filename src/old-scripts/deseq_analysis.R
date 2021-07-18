library(tximport)
library(data.table)
library(DESeq2)
library(ggplot2)
library(dplyr)

gene2tx <- fread("data/current_assembly/output/trinity/Trinity.fasta.gene_trans_map", header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files (look in output/salmon for files named quant.sf)
  quant_files <- list.files(path="output/2020_assembly_output/mh_salmon/", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
##assign names to quant files from folder name
  names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)
##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
##for methods that only provide transcript level estimates e.g. salmon)
  txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
##Make a data frame describing samples (here the only info is name and sample which happen to be the same)
  sample_data <- data.frame(row.names = names(quant_files), sample=names(quant_files))

##design is where to chuck in more info about how you want to analyse data(?)
  dds <- DESeqDataSetFromTximport(txi, colData = sample_data, design = ~sample)
  saveRDS(dds, "output/2019_assembly_ovaries_output/deseq2/dds.rds")
  dds <- readRDS("output/2019_assembly_ovaries_output/deseq2/dds.rds")
##get table of counts for each gene - what samples show expression of lbfv genes?
  counts_table <- data.table((data.frame(counts(dds))), keep.rownames = TRUE)
  counts_matrix <- data.table(data.frame(counts(dds)), keep.rownames = TRUE)
  counts_colSums <- data.frame(colSums(counts_matrix, na.rm=TRUE))
  
annots <- fread("data/current_assembly/output/trinotate/sorted/best_annot_per_gene.csv")
  
##read in list of transcripts that matched to viral scaffolds
  transcripts_viral_scaffolds <- fread("data/2019_with_ovaries/transcripts_annots_sorted.csv", na.strings = "")
  transcripts_viral_scaffolds$gene_id <- tstrsplit(transcripts_viral_scaffolds$transcript_id, "_i", keep=c(1))
  recip_blast_viral <- fread("data/2019_with_ovaries/recip_blastx_viral_transcripts.csv")
  recip_blast_viral$gene_id <- tstrsplit(recip_blast_viral$transcript_id, "_i", keep=c(1))
  viral_transcripts_annots <- merge(transcripts_viral_scaffolds, recip_blast_viral, by="gene_id", all=TRUE)
  viral_transcripts_annots <- viral_transcripts_annots[,c(1,3,5)]
  setnames(viral_transcripts_annots, old=c("annotation.x", "annotation.y"), new=c("genome_hit_annot", "recip_transcript_blast_annot"))
  trinotate_viral_annots <- fread("data/2019_with_ovaries/gene_id_vs_viral_annot.csv")
  viral_transcripts_annots <- merge(viral_transcripts_annots, trinotate_viral_annots, by.x="gene_id", by.y="#gene_id", all=TRUE)
  
  ##merge with counts table
  viral_scaffold_transcript_counts <- merge(viral_transcripts_annots, counts_table, by.x="gene_id", by.y="rn", all.x=TRUE)
  fwrite(viral_scaffold_transcript_counts, "output/2019_assembly_ovaries_output/deseq2/deseq_counts_viral_scaffold_genes.csv")
    
  ##what about transcripts with viral trinotate hits?
  
########################################  
##no replicates, can't do this anymore##
######################################## 
  dds <- DESeq(dds)

###sting vs head
  ##extract a results table using contrast to specify your comparison (look up for more detail)
    stingvheadres <- results(dds, contrast = (c("sample", "sting", "head")))
  ##make nice data table out of results
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

  
###abdo vs head
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


###sting vs abdo
    res <- results(dds, contrast = (c("sample", "sting", "abdo")))
    results_table <- data.table(data.frame(res), keep.rownames = TRUE)
    setkey(results_table, log2FoldChange)
    results_table[!is.na(padj)]
    fwrite(results_table, "output/deseq2/sting_vs_abdo_results.txt")

    
##other playing around
  
##make csv of sting quant results
sting <- tximport("output/salmon/sting_quant/quant.sf", type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
stingquant<- data.frame(sting)
write.csv(stingquant, file = "output/new_transcriptome_deseq2/stingquant.csv")

##Plot genes of intrest after blast search - change gene= to a gene of interest
d <- plotCounts(dds, gene="TRINITY_DN21692_c3_g3", intgroup="sample", returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
  
###????? 
##To subset annotation report into significant genes
annotation.report  <- fread("output/annotation_report.csv")
sig.annotations <- annotation.report[gene_id %in% stingvhead_resSig$V1]
fwrite(sig.annotations, "output/new_transcriptome_deseq2/stingvhead_sigannotations.csv")
