library(data.table)
library(tidyverse)
library(dplyr)
library(pheatmap)
library(viridis)

trinotate <- fread("data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv")

odorant_bx <- subset(trinotate, grepl("odorant", sprot_Top_BLASTX_hit, ignore.case=TRUE))
odorant_bp <- subset(trinotate, grepl("odorant", sprot_Top_BLASTP_hit, ignore.case=TRUE))
odorant <- full_join(odorant_bx, odorant_bp)

##ORs
odorant_receptor_bx <- subset(odorant, grepl("receptor", sprot_Top_BLASTX_hit, ignore.case=TRUE))
odorant_receptor_bp <- subset(odorant, grepl("receptor", sprot_Top_BLASTP_hit, ignore.case=TRUE))
odorant_receptor <- full_join(odorant_receptor_bx, odorant_receptor_bp)
odorant_receptor$Label <- paste("OR")
fwrite(odorant_receptor, "output/deseq2/tissue_itWT_LRT/head/chemosensory/ORs.csv")


##OBPs
odorant_binding_bx <- subset(odorant, grepl("binding", sprot_Top_BLASTX_hit, ignore.case=TRUE))
odorant_binding_bp <- subset(odorant, grepl("binding", sprot_Top_BLASTP_hit, ignore.case=TRUE))
odorant_binding <- full_join(odorant_binding_bx, odorant_binding_bp)
odorant_binding$Label <- paste("OBP")

##GRs
gustatory_receptor_bx <- subset(trinotate, grepl("gustatory receptor|gustatory and pheromone receptor", sprot_Top_BLASTX_hit, ignore.case=TRUE))
gustatory_receptor_bp <- subset(trinotate, grepl("gustatory receptor|gustatory and pheromone receptor", sprot_Top_BLASTP_hit, ignore.case=TRUE))
gustatory_receptor <- full_join(gustatory_receptor_bx, gustatory_receptor_bp)
gustatory_receptor$Label <- paste("GR")

##IRs
ionotropic_receptor_bx <- subset(trinotate, grepl("ionotropic receptor", sprot_Top_BLASTX_hit, ignore.case=TRUE))
ionotropic_receptor_bp <- subset(trinotate, grepl("ionotropic receptor", sprot_Top_BLASTP_hit, ignore.case=TRUE))
ionotropic_receptor <- full_join(ionotropic_receptor_bx, ionotropic_receptor_bp)
ionotropic_receptor$Label <- paste("IR")

##SNMPs
snmp_bx <- subset(trinotate, grepl("sensory neuron membrane protein", sprot_Top_BLASTX_hit, ignore.case=TRUE))
snmp_bp <- subset(trinotate, grepl("sensory neuron membrane protein", sprot_Top_BLASTP_hit, ignore.case=TRUE))
snmp <- full_join(snmp_bx, snmp_bp)
snmp$Label <- paste("SNMP")

##all chemosensory genes
chemosensory <- full_join(odorant_binding, full_join(odorant_receptor, full_join(gustatory_receptor, full_join(ionotropic_receptor, snmp))))
fwrite(chemosensory, "output/deseq2/tissue_itWT_LRT/head/chemosensory/all_chemosensory.csv")

##merge with salmon mean TPM?
tpm_means <- fread("output/deseq2/salmon_TPM_means.csv")
chemosensory_tpm <- merge(chemosensory, tpm_means, by.x="#gene_id", by.y="rn")

lrt_res <- fread("output/deseq2/tissue_LRT/sig_degs.csv")
chemosensory_degs <- merge(chemosensory, lrt_res, by.x="#gene_id", by.y="rn")
fwrite(chemosensory_degs, "output/deseq2/tissue_itWT_LRT/head/chemosensory/chemosensory_lrt_degs.csv")

id_to_group_dt <- chemosensory_degs[,c(1,18)]
id_to_group <- id_to_group_dt %>% remove_rownames %>% column_to_rownames(var="#gene_id")

## heatmap ##
mh_dds <- readRDS("output/deseq2/tissue_LRT/all_mh_tissue_LRT.rds")
##vst transform
mh_vst <- varianceStabilizingTransformation(mh_dds, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)
##subset for DEGs
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% chemosensory_degs$`#gene_id`)
##reorder rows based on group
plot_group <- merge(mh_vst_degs, id_to_group_dt, by.x="rn", by.y="#gene_id")
ordered <- plot_group[order(Label, -Mh_head1),]
mh_vst_degs_plot <- ordered[,-c(20)]

##turn first row back to row name
mh_vst_degs_plot <- mh_vst_degs_plot %>% remove_rownames %>% column_to_rownames(var="rn")
##reorder tissues for plot
mh_vst_degs_plot <- mh_vst_degs_plot[,c(4,5,6,16,10,11,12,1,2,3,15,7,8,9,17,18,13,14)]


##get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mh_dds)[,c("Tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
##for plot label
sample_to_tissue <- as.data.frame(colData(mh_dds)[,c("Tissue", "sample_name")])
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")

tissue_colours <- list(Tissue = c(Head="#440154FF", Thorax="#3B528BFF", Abdomen="#21908CFF", Ovaries="#5DC863FF", Venom="#FDE725FF"),
                       Label = c(GR="#3B0F70FF", IR="#8C2981FF",  OBP="#DE4968FF", OR="#FE9F6DFF", SNMP="#FCFDBFFF"))

##plot
pheatmap(mh_vst_degs_plot, cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, annotation_row=id_to_group,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))

