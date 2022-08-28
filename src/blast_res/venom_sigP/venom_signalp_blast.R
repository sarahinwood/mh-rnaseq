library(data.table)
library(dplyr)
library(tidyr)

signalp_blast <- fread("output/04_blast/venom_trinotate_signalp/venom_signalp_nr_blastx.outfmt6")

setnames(signalp_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"),
         new=c("transcript_id", "hit", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(signalp_blast, transcript_id, evalue, -bit_score)
signalp_blast$`#gene_id` <- tstrsplit(signalp_blast$transcript_id, "_i", keep=c(1))

##filter out hypothetical/unchar
no_hypo <- subset(signalp_blast, !grepl("uncharacter|hypothetical|unnamed|unknown|GSCOCG", annotation, ignore.case = TRUE))

##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues <- no_hypo[,.SD[which.min(evalue)], by=transcript_id]
##write table
fwrite(min_evalues, "output/04_blast/venom_trinotate_signalp/venom_trinotate_signalp_besthits.csv")

venom_degs <- fread("output/03_deseq2/tissue_itWT_LRT/Venom/Venom_sp_LRT_all_annots.csv")
all_res <- merge(venom_degs, min_evalues, by="transcript_id", all=T)
fwrite(all_res, "output/03_deseq2/tissue_itWT_LRT/Venom/Venom_sp_LRT_sigp_blast.csv")
