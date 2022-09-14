#!/usr/bin/env Rscript
#############
# LIBRARIES #
#############

library(data.table)
library(tidyverse)

########
# MAIN #
########

blast_res <- fread("output/03_deseq2/tissue_itWT_LRT/meiosis_sp/dmc1_rec8/dmc1_rec8_nr_blastx.outfmt6")
trinotate <- fread('data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv', na.strings="")

setnames(blast_res, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"),
         new=c("transcript_id", "hit", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))


##filter out hypothetical/unchar
unann_blast_no_hypo <- subset(blast_res, !grepl("uncharacter|hypothetical|unnamed|unknown|GSCOCG", annotation, ignore.case = TRUE))

##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(unann_blast_no_hypo, transcript_id, evalue, -bit_score)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues <- unann_blast_no_hypo[,.SD[which.min(evalue)], by=transcript_id]
dmc1_rec8 <- subset(min_evalues, grepl("DMC1|REC8", annotation, ignore.case = TRUE))
dmc1_rec8_trinotate <- merge(dmc1_rec8, trinotate, by="transcript_id")

fwrite(min_evalues, "output/03_deseq2/tissue_itWT_LRT/meiosis_sp/dmc1_rec8/dmc1_rec8_hits.csv")