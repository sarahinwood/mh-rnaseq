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

#############
## GLOBALS ##
#############

venom_file <- snakemake@input[["venom_file"]]

########
# MAIN #
########

venom_degs <- fread(venom_file, na.strings = "")
# subset for signalP
signalp <- subset(venom_degs, SignalP=="Y")
sigp_no_nr <- subset(signalp, is.na(annotation))

fwrite(list(sigp_no_nr$transcript_id), snakemake@output[["deg_ids"]])

# write log
sessionInfo()