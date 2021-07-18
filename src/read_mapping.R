library(data.table)
library(DESeq2)

##mh reads mapped
mh_dds <- readRDS("output/deseq2/mh_dds.rds")
counts_table_mh <- (data.table(counts(mh_dds)))
counts_colSums_mh <- setDT(data.frame(colSums(counts_table_mh, na.rm=TRUE)), keep.rownames=TRUE)

setnames(counts_colSums_mh, old=c("rn", "colSums.counts_table_mh..na.rm...TRUE."), new=c("Sample_name", "readpairs_mapped_mh"))

##both reads mapped
bbduk_reads_out <- fread("output/bbduk_reads_out.csv")
full_read_mapping <- merge(counts_colSums_mh, bbduk_reads_out, by="Sample_name")
full_read_mapping$bbduk_halved <- (full_read_mapping$bbduk_reads_out)/2
full_read_mapping$`%_ofall_mh` <- (full_read_mapping$readpairs_mapped_mh/full_read_mapping$bbduk_halved)*100
##slight differences from salmon - due to salmon filtering out some extra reads first
fwrite(full_read_mapping, "output/deseq2/read_mapping.csv")
