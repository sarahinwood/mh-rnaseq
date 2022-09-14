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

trinotate_file <- snakemake@input[["trinotate_file"]]

########
# MAIN #
########

##don't miss any by using length-filtered
trinotate_report <- fread(trinotate_file)

meiosis_bx <- subset(trinotate_report, grepl("CORTEX|SPO11|HOP2|MND1|DMC1|MSH4|MSH5|REC8", sprot_Top_BLASTX_hit, ignore.case=TRUE))
meiosis_bp <- subset(trinotate_report, grepl("CORTEX|SPO11|HOP2|MND1|DMC1|MSH4|MSH5|REC8", sprot_Top_BLASTP_hit, ignore.case=TRUE))
meiosis_pfam <- subset(trinotate_report, grepl("CORTEX|SPO11|HOP2|MND1|DMC1|MSH4|MSH5|REC8|MUTS", Pfam, ignore.case=TRUE))
# filter pfam hits
meiosis_pfam_fil <- subset(meiosis_pfam, !grepl("MSH2_HUMAN|MSH6_DROME|MUTS_THEMA|RAD21_BOVIN", sprot_Top_BLASTX_hit, ignore.case=TRUE))
#TRINITY_DN22299_c0_g1 DMC1, TRINITY_DN85877_c0_g1 REC8
dmc1_rec8_from_blast <- subset(trinotate_report, grepl("TRINITY_DN22299_c0_g1|TRINITY_DN85877_c0_g1", `#gene_id`, ignore.case=TRUE))
trinotate_meiosis_sp <- full_join(meiosis_bx, full_join(meiosis_bp, full_join(meiosis_pfam_fil, dmc1_rec8_from_blast)))

# label gene name
trinotate_meiosis_sp$bx_gene <- tstrsplit(trinotate_meiosis_sp$sprot_Top_BLASTX_hit, "_", keep=c(1))
trinotate_meiosis_sp$bp_gene <- tstrsplit(trinotate_meiosis_sp$sprot_Top_BLASTP_hit, "_", keep=c(1))
trinotate_meiosis_sp$pfam_gene <- tstrsplit(trinotate_meiosis_sp$Pfam, "^", keep=c(3), fixed=T)
# add gene label
trinotate_meiosis_sp$meiosis_gene <- ifelse(trinotate_meiosis_sp$`#gene_id` %in% meiosis_bx$`#gene_id`, paste(trinotate_meiosis_sp$bx_gene),
                                            ifelse(trinotate_meiosis_sp$`#gene_id` %in% meiosis_bp$`#gene_id`, paste(trinotate_meiosis_sp$bp_gene),
                                                   ifelse(trinotate_meiosis_sp$`#gene_id` %in% meiosis_pfam$`#gene_id`, paste(trinotate_meiosis_sp$pfam_gene),
                                                          ifelse(trinotate_meiosis_sp$`#gene_id`=="TRINITY_DN22299_c0_g1", paste("DMC1"), paste("REC8")))))

fwrite(trinotate_meiosis_sp, snakemake@output[["meiosis_sp_genes"]])
fwrite(list(trinotate_meiosis_sp$transcript_id), snakemake@output[["meiosis_ids"]])


# write log
sessionInfo()