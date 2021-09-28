library(data.table)
library(dplyr)

crawford_blast_res <- fread("output/blast/crawford_venom/blastn.outfmt6")
trinotate <- fread("data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv", na.strings = ".")

setnames(crawford_blast_res, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"),
         new=c("Crawford_seq", "Trinity_ID", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(crawford_blast_res, Crawford_seq, evalue, -bit_score)
crawford_blast_res$`#gene_id` <- tstrsplit(crawford_blast_res$Trinity_ID, "_i", keep=c(1))
crawford_blast_res$short_gene_id <- tstrsplit(crawford_blast_res$`#gene_id`, "TRINITY_DN", keep=c(2))
##table of all transcripts with crawford hits
fwrite(crawford_blast_res, "output/blast/crawford_venom/crawford_all_hits.csv")

##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues <- crawford_blast_res[,.SD[which.min(evalue)], by=Crawford_seq]
##merge with trinotate results
min_evalues_annots <- merge(min_evalues, trinotate, by.x="Trinity_ID", by.y="transcript_id")
##write table
fwrite(min_evalues_annots, "output/blast/crawford_venom/crawford_best_hits.csv")

#overlap with venom DEGs
venom_degs <- fread("output/deseq2/tissue_itWT_LRT/venom/venom_sp_LRT_annots.csv")
venom_intersect <- intersect(crawford_blast_res$Trinity_ID, venom_degs$transcript_id)
venom_intersect_crawford_blast <- merge(venom_intersect_degs, crawford_blast_res, by.x="transcript_id", by.y="Trinity_ID")
venom_intersect_crawford_blastdt <- venom_intersect_crawford_blast[,c(1,2,3,4,8,9,14,25,33,35,36,37,45)]
fwrite(venom_intersect_crawford_blastdt, "output/blast/crawford_venom/venom_transcripts_annots_degs.csv")

############################
## plot venom gene counts ##
############################
##plot venom genes not DE
mh_dds_lrt <- readRDS("output/deseq2/stage_WT/mh_stage_WT.rds")
##get gene counts
counts_table <- data.table(counts(mh_dds_lrt, normalized=TRUE), keep.rownames = TRUE)
annot_counts <- filter(counts_table, rn %in% crawford_blast_res$`#gene_id`)
##melt for plotting
plot_annots_counts <- annot_counts %>% gather(colnames(annot_counts)[2:22], key="sample_name", value="normalized_counts")
##sample group information
sample_table <- fread("data/sample_table.csv")
sample_table$group <- paste(sample_table$tissue)
name_vs_group <- sample_table[,c(1,2)]
plotting_counts <- inner_join(plot_annots_counts, name_vs_group)
group_order <- c("Head", "Thorax", "Abdomen", "Ovaries", "Venom", "Pupa")
plotting_counts$group <- factor(plotting_counts$tissue, levels=group_order)
##merge with crawford gene annot
crawford_blast_res$label <- tstrsplit(crawford_blast_res$Crawford_seq, "hyperodae_", keep=c(2))
crawford_blast_res$label <- tstrsplit(crawford_blast_res$label, "_mRNA", keep=c(1))
gene_to_crawford <- crawford_blast_res[,c(14,15,16)]
plotting_counts_crawford_hits <- merge(plotting_counts, gene_to_crawford, by.x="rn", by.y="#gene_id")
plotting_counts_crawford_hits$plot_label <- paste(plotting_counts_crawford_hits$short_gene_id, plotting_counts_crawford_hits$label, sep="_")
##plot all annot DEGs using ggplot2
ggplot(plotting_counts_crawford_hits) +
  geom_point(aes(x = group, y = normalized_counts, colour=group)) +
  labs(colour="Tissue", y="Normalized counts", x="")+
  scale_colour_viridis(discrete=TRUE)+
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~plot_label, scales="free")

##for best hits##
counts_table <- data.table(counts(mh_dds_lrt, normalized=TRUE), keep.rownames = TRUE)
best_annot_counts <- filter(counts_table, rn %in% min_evalues_annots$`#gene_id.x`)
##melt for plotting
plot_annots_counts <- annot_counts %>% gather(colnames(annot_counts)[2:22], key="sample_name", value="normalized_counts")
##sample group information
sample_table <- fread("data/sample_table.csv")
sample_table$group <- paste(sample_table$tissue)
name_vs_group <- sample_table[,c(1,2)]
plotting_counts <- inner_join(plot_annots_counts, name_vs_group)
group_order <- c("Head", "Thorax", "Abdomen", "Ovaries", "Venom", "Pupa")
plotting_counts$group <- factor(plotting_counts$tissue, levels=group_order)
##merge with crawford gene annot
min_evalues_annots$label <- tstrsplit(min_evalues_annots$Crawford_seq, "hyperodae_", keep=c(2))
min_evalues_annots$label <- tstrsplit(min_evalues_annots$label, "_mRNA", keep=c(1))
gene_to_crawford <- min_evalues_annots[,c(14,32)]
plotting_counts_crawford_hits <- merge(plotting_counts, gene_to_crawford, by.x="rn", by.y="#gene_id.x")
##plot all annot DEGs using ggplot2
ggplot(plotting_counts_crawford_hits) +
  geom_point(aes(x = group, y = normalized_counts, colour=group)) +
  labs(colour="Tissue", y="Normalized counts", x="")+
  scale_colour_viridis(discrete=TRUE)+
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~label, scales="free")

