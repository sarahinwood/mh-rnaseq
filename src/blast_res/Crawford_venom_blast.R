library(data.table)
library(dplyr)

crawford_blast_res <- fread("output/blast/crawford_venom/blastn.outfmt6")
trinotate <- fread("data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings = ".")

setnames(crawford_blast_res, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"),
         new=c("Crawford_seq", "Trinity_ID", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(crawford_blast_res, Crawford_seq, evalue, -bit_score)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues <- crawford_blast_res[,.SD[which.min(evalue)], by=Crawford_seq]

##merge with trinotate results
min_evalues_annots <- merge(min_evalues, trinotate, by.x="Trinity_ID", by.y="transcript_id")
##write table
fwrite(min_evalues_annots, "output/blast/crawford_venom/crawford_best_hits.csv")
