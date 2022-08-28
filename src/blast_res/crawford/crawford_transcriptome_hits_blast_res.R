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

library(data.table)
library(tidyverse)

#############
## GLOBALS ##
#############

res_file <- snakemake@input[["res_file"]]

########
# MAIN #
########

blast_res <- fread(res_file)
setnames(blast_res, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"),
         new=c("transcript_id", "hit", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
blast_res <- subset(blast_res, !grepl("uncharacter|hypothetical", annotation))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(blast_res, transcript_id, evalue, -bit_score)
##table of all transcripts with crawford hits
fwrite(blast_res, snakemake@output[["all_res"]])

#remove hyperodae results
blast_res_no_hyp <- subset(blast_res, !grepl("hyperodae", annotation))
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues <- blast_res_no_hyp[,.SD[which.min(evalue)], by=transcript_id]

##write table
fwrite(min_evalues, snakemake@output[["best_hits"]])


# write log
sessionInfo()

