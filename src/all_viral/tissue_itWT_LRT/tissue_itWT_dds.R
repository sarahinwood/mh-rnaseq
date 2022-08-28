library(tximport)
library(data.table)
library(DESeq2)

mh_dds <- readRDS("output/03_deseq2/mh_dds.rds")

viral_genes <- fread("data/mh-all-transcriptome/output/recip_blast/viral_nr_blastx/best_viral_hits_plot.csv")
viral_genes$gene_id <- tstrsplit(viral_genes$transcript_id, "_i", keep=c(1))

##select only adult tissue samples and only viral genes
mh_dds_viral <- mh_dds[viral_genes$gene_id,]
mh_dds_tissue <- mh_dds_viral[,mh_dds_viral$stage=="adult"]
##re-level tissue factor
mh_dds_tissue$Tissue <- factor(mh_dds_tissue$tissue, levels=c("Head", "Thorax", "Abdomen", "Ovaries", "Venom"))
mh_dds_tissue$Rep <- factor(mh_dds_tissue$rep)
mh_dds_tissue$Flowcell <- factor(mh_dds_tissue$flowcell)

design(mh_dds_tissue) <- ~Flowcell+Rep+Tissue
mh_dds_tissue <- DESeq(mh_dds_tissue)
saveRDS(mh_dds_tissue, "output/03_deseq2/viral_LRT/viral_itWT.rds")
