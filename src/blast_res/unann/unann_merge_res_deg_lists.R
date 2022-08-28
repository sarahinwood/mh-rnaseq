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

###########
# GLOBALS #
###########

blast_res_file <- snakemake@input[["blast_res_file"]]
deg_file <- snakemake@input[["deg_file"]]

########
# MAIN #
########

blast_res <- fread(blast_res_file)
degs <- fread(deg_file, na.strings="")
# merge
degs_all_annots <- merge(degs, blast_res, all.x=T)

fwrite(degs_all_annots, snakemake@output[["degs_all_annots"]])

# write log
sessionInfo()
