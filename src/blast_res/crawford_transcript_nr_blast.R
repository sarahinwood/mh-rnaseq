library(data.table)
library(dplyr)

blast_res <- fread("output/blast/crawford_venom/crawford_transcripts.outfmt6")
setnames(blast_res, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"),
         new=c("Transcript_ID", "hit", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
blast_res_not_hyp <- subset(blast_res, !grepl("hyperodae", annotation))
blast_res_not_hyp <- subset(blast_res_not_hyp, !grepl("uncharacter", annotation))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(blast_res_not_hyp, Transcript_ID, evalue, -bit_score)
##table of all transcripts with crawford hits
fwrite(blast_res_not_hyp, "output/blast/crawford_venom/transcripts_all_nrblast.csv")

##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues <- blast_res_not_hyp[,.SD[which.min(evalue)], by=Transcript_ID]

##write table
fwrite(min_evalues, "output/blast/crawford_venom/transcripts_best_nrblast.csv")
