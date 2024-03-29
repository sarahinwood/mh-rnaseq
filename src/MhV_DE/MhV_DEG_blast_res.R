library(data.table)
library(dplyr)

unann_blast <- fread("output/deseq2/MhV_LRT/MhV_degs.outfmt6")

##viral_contig_blast res
setnames(unann_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"),
         new=c("transcript_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))

##filter out hypothetical/unchar?
unann_blast_no_pieris <- subset(unann_blast, !grepl("Pieris macdunnoughi|CINCED", annotation, ignore.case = TRUE))

##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(unann_blast_no_pieris, transcript_id, evalue, -bit_score)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues <- unann_blast_no_pieris[,.SD[which.min(evalue)], by=transcript_id]

fwrite(min_evalues, "output/deseq2/MhV_LRT/min_evalues.csv")
