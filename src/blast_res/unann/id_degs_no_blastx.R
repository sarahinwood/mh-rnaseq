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

########
# MAIN #
########

tissue_files <- list.files("output/03_deseq2/tissue_itWT", "_sp_annots.csv", recursive=T, full.names = T)
stage_file <- "output/03_deseq2/stage_WT/sig_w_blast_annots.csv"
deg_files <- c(tissue_files, stage_file)

all_degs <- lapply(deg_files, read_csv) %>% bind_rows
no_blastx <- subset(all_degs, is.na(sprot_Top_BLASTX_hit))
degs_no_blastx <- unique(no_blastx$rn)

fwrite(list(degs_no_blastx), snakemake@output[["deg_ids_no_blastx"]])

# write log
sessionInfo()
