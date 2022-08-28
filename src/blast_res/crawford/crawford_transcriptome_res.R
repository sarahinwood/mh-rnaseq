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
trinotate_file <- snakemake@input[["trinotate_file"]]
venom_degs_file <- snakemake@input[["venom_degs_file"]]
transcript_lengths_file <- snakemake@input[["transcript_lengths_file"]]
crawford_lengths_file <- snakemake@input[["crawford_lengths_file"]]
salmon_tpm_file <- snakemake@input[["salmon_tpm_file"]]

########
# MAIN #
########

crawford_blast_res <- fread(res_file, select=c(1:12))
trinotate <- fread(trinotate_file, na.strings = ".")
transcript_lengths <- fread(transcript_lengths_file, select=c(1,3), col.names=c("transcript_id", "transcript_length"))
crawford_lengths <- fread(crawford_lengths_file)
salmon_tpm <- fread(salmon_tpm_file)

## set up salmon mean tpm table
rest_mean_tpm <- data.frame(ID=salmon_tpm[,1], rest_means=rowMeans(salmon_tpm[,2:16, 19:22]))
venom_mean_tpm <- data.frame(ID=salmon_tpm[,1], venom_means=rowMeans(salmon_tpm[,17:18]))
setorder(venom_mean_tpm, -venom_means)
# venom tpm rank
venom_mean_tpm$venom_rank <- seq.int(nrow(venom_mean_tpm))
all_tpm <- merge(rest_mean_tpm, venom_mean_tpm, by="rn")

## blast res
setnames(crawford_blast_res, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"),
         new=c("crawford_gene", "transcript_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(crawford_blast_res, crawford_gene, evalue, -bit_score)
crawford_blast_res$`#gene_id` <- tstrsplit(crawford_blast_res$transcript_id, "_i", keep=c(1))
crawford_blast_res$short_gene_id <- tstrsplit(crawford_blast_res$`#gene_id`, "TRINITY_DN", keep=c(2))
# merge with crawford & transcript lengths
crawford_blast_res_length <- merge(crawford_blast_res, crawford_lengths, by="crawford_gene")
crawford_blast_res_length <- merge(crawford_blast_res_length, transcript_lengths, by="transcript_id", all.x=T)
##table of all transcripts with crawford hits
fwrite(crawford_blast_res_length, snakemake@output[["all_res"]])

# remove chimeric VG5 transcript hit - TRINITY_DN1215_c1_g1_i10
crawford_blast_res_length_fil <- subset(crawford_blast_res_length, transcript_id!="TRINITY_DN1215_c1_g1_i10")
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie - protein 5 has CO hit again - chimeric? has another GILT that is DE
min_evalues <- crawford_blast_res_length_fil[,.SD[which.min(evalue)], by=crawford_gene]
##merge with trinotate results
min_evalues_annots <- merge(min_evalues, trinotate, by="transcript_id")
##write table
fwrite(min_evalues_annots, snakemake@output[["best_res"]])
fwrite(list(min_evalues_annots$transcript_id), snakemake@output[["best_hit_transcript_ids"]])

#overlap with venom DEGs
venom_degs <- fread(venom_degs_file)
venom_intersect <- intersect(crawford_blast_res$transcript_id, venom_degs$transcript_id)
venom_intersect_crawford_blast <- merge(venom_degs, crawford_blast_res, by="transcript_id", all.y=T)
venom_intersect_crawford_blast_degs <- subset(venom_intersect_crawford_blast, !is.na(padj))
venom_intersect_crawford_blast_degs$venom_protein <- tstrsplit(venom_intersect_crawford_blast_degs$crawford_gene, "hyperodae_", keep=c(2))
venom_intersect_crawford_blast_degs$venom_protein <- tstrsplit(venom_intersect_crawford_blast_degs$venom_protein, "_mRNA", keep=c(1))
degs_tpms <- merge(venom_intersect_crawford_blast_degs, all_tpm, by="rn")
degs_tpms_dt <- degs_tpms[,c(37,1:4,8,9,14,15,24,25,33,34, 38:40)]
fwrite(degs_tpms_dt, snakemake@output[["crawford_venom_degs"]])

# write log
sessionInfo()


