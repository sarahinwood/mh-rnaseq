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

library(tidyverse)
library(data.table)


###########
# GLOBALS #
###########

blast_res_file <- snakemake@input[["blast_res_file"]]

########
# MAIN #
########

unann_blast <- fread(blast_res_file)

##viral_contig_blast res
setnames(unann_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"),
         new=c("transcript_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))

##filter out hypothetical/unchar
unann_blast_no_hypo <- subset(unann_blast, !grepl("uncharacter|hypothetical|unnamed|unknown|GSCOCG", annotation, ignore.case = TRUE))

##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(unann_blast_no_hypo, transcript_id, evalue, -bit_score)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues <- unann_blast_no_hypo[,.SD[which.min(evalue)], by=transcript_id]

fwrite(min_evalues, snakemake@output[["blast_best_res"]])

# write log
sessionInfo()
