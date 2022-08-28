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
library(ggplot2)
library(viridis)
library(cowplot)

###########
# GLOBALS #
###########

##DEG_list             
pfam_file <- snakemake@input[["pfam_file"]]
go_file <- snakemake@input[["go_file"]]

########
# MAIN #
########

pfam_enrich <- fread(pfam_file)
pfam_enrich$mod <- ifelse(pfam_enrich$ID=="PF00732.19", "N",
                              ifelse(pfam_enrich$ID=="PF05199.13", "C",
                                     ifelse(pfam_enrich$ID=="PF13895.6", "2", 
                                            ifelse(pfam_enrich$ID=="PF13927.6", "3", ""))))
pfam_enrich$Description <- paste(pfam_enrich$Description, pfam_enrich$mod, sep=" ")
pfam_enrich$Description <- factor(pfam_enrich$Description, levels=pfam_enrich$Description[order(pfam_enrich$GeneRatio, pfam_enrich$Description, decreasing=F)])
go_enrich <- fread(go_file)
go_enrich$Description <- factor(go_enrich$Description, levels=go_enrich$Description[order(go_enrich$pathway_kind, go_enrich$GeneRatio, go_enrich$Description, decreasing=T)])

pfam <- ggplot(pfam_enrich, aes(x=Description, y=GeneRatio)) +
  geom_col(aes(fill="#440154FF"))+
  labs(x="Pfam domain", y="GeneRatio", size="Leading\nedge size") +
  coord_flip() +
  scale_y_continuous(breaks=c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14), limits=c(0, 0.14))+
  scale_fill_viridis(discrete=TRUE)+
  theme_bw()+
  theme(legend.position = "none")


go <- ggplot(go_enrich, aes(Description, GeneRatio)) +
  geom_col(aes(fill=pathway_kind))+
  labs(x="Gene ontology terms", y="GeneRatio",
       fill="GO domain", size="Leading\nedge size") +
  coord_flip() +
  scale_y_continuous(breaks=c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14), limits=c(0, 0.14))+
  scale_fill_viridis(discrete=TRUE, begin=0.5)+
  theme_bw()+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

go_length <- length(go_enrich$Description)
pfam_length <- length(pfam_enrich$Description)
go_plot_height <-go_length/pfam_length

pdf(snakemake@output[["enrich_plot"]], width=7, height=5) # default 7x7
plot_grid(go, pfam, ncol=1, align="v",
          rel_heights=c(go_plot_height, 1))
dev.off()

# write log
sessionInfo()
