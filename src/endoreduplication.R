library(data.table)

trinotate_report <- fread("data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv")
ovary_DEGs <- fread("output/deseq2/tissue_itWT_LRT/ovary/ovary_sp_LRT_annots.csv")

##4 BlastX GO terms associated with endoreduplication (no Pfam GO are, BlastP has less)
###GO:0042023 DNA endoreduplication
###GO:0032877 positive regulation of DNA endoreduplication
###GO:0032876 negative regulation of DNA endoreduplication
###GO:0032875 regulation of DNA endoreduplication

endo_GObx <- subset(trinotate_report, grepl("endoreduplication", gene_ontology_BLASTX, ignore.case=TRUE))
endo_GObp <- subset(trinotate_report, grepl("endoreduplication", gene_ontology_BLASTP, ignore.case=TRUE))
endo_GO <-full_join(endo_GObx, endo_GObp)

##endo genes DE in LRT test
lrt_degs <- fread("output/deseq2/tissue_itWT_LRT/ovary/sig_ovary_lrt.csv")
endo_degs <- subset(lrt_degs, rn %in% endo_GO$`#gene_id`)
endo_degs_annots <- merge(endo_degs, trinotate_report, by.x="rn", by.y="#gene_id")
endo_degs_annots$gene <- tstrsplit(endo_degs_annots$sprot_Top_BLASTX_hit, "_", keep=c(1))
##* for ovary specific also
endo_degs_annots$DE <- ifelse(endo_degs_annots$rn %in% ovary_DEGs$rn, "*", "")
endo_degs_annots$gene_de <- paste(endo_degs_annots$gene, endo_degs_annots$DE, sep="")
id_to_gene <- endo_degs_annots[,c(1,26)]

##add GO label for negative, positive, endoreduplication
negative <- subset(trinotate_report, grepl("GO:0032876", gene_ontology_BLASTX, ignore.case=TRUE))
positive <- subset(trinotate_report, grepl("GO:0032877", gene_ontology_BLASTX, ignore.case=TRUE))
endoreduplication <- subset(trinotate_report, grepl("GO:0042023", gene_ontology_BLASTX, ignore.case=TRUE))

##GO label
gene_id_GO_label <- id_to_gene
gene_id_GO_label$plot_gene_id <- paste(gene_id_GO_label$gene_de, gene_id_GO_label$rn, sep="_")
gene_id_GO_label$GO <- ifelse(gene_id_GO_label$rn %in% negative$`#gene_id`, "Negative regulation",
                              ifelse(gene_id_GO_label$rn %in% positive$`#gene_id`, "Positive regulation", "Endoreduplication"))
gene_id_GO_label <- gene_id_GO_label[,c(3,4)]
gene_id_GO_label <- gene_id_GO_label %>% remove_rownames %>% column_to_rownames(var="plot_gene_id")



##endo DEGs in ovary specific
ovary_endo <- subset(ovary_DEGs, rn %in% endo_GO$`#gene_id`)

#############
## heatmap ##
#############

mh_dds_lrt <- readRDS("output/deseq2/tissue_itWT_LRT/ovary/ovary_LRT_dds.rds")
##vst transform
mh_vst <- varianceStabilizingTransformation(mh_dds_lrt, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)

##subset for genes of interest
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% id_to_gene$rn)
mh_vst_degs_label <- merge(mh_vst_degs, id_to_gene, by="rn")
mh_vst_degs_label$rn <- paste(mh_vst_degs_label$gene_de, mh_vst_degs_label$rn, sep="_")
mh_vst_degs_label <- mh_vst_degs_label %>% remove_rownames %>% column_to_rownames(var="rn")
##reorder for plot - old samples are off to side for now - need to change this !!!!!!!!!!!!!!!!!!!!!
mh_vst_degs_plot <- mh_vst_degs_label[,c(7,8,9,17,18,4,5,6,16,10,11,12,1,2,3,15,13,14)]

##get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
##for plot label
sample_to_tissue <- as.data.frame(colData(mh_dds_lrt)[,c("Tissue", "sample_name")])
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")

tissue_colours <- list(Tissue = c(Head="#440154FF", Thorax="#3B528BFF", Abdomen="#21908CFF", Ovaries="#5DC863FF", Venom="#FDE725FF"),
                       GO = c("Negative regulation"="#51127CFF", "Positive regulation"="#B63679FF", "Endoreduplication"="#FB8861FF"))
##plot
##not clustered by sample
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=TRUE, annotation_row=gene_id_GO_label,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))

#open PDF in inkscape and change gene labels to remove trinity IDs, bold sig DEGs

