library(data.table)
library(DESeq2)

#highly expressed venom genes
salmon_tpm <- fread("output/03_deseq2/salmon_mean_TPM.csv")
venom_500 <- subset(salmon_tpm, venom_rank<=500)
# remove those DE vs all tissues
venom_degs <- fread("output/03_deseq2/tissue_itWT_LRT/Venom/Venom_sp_LRT_sigp_blast.csv")
venom_500_notDE <- subset(venom_500, !(venom_500$rn %in% venom_degs$rn))
# signal peptides?
trinotate <- fread("data/mh-all-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv")
venom_500_notDE_sigp <- merge(venom_500_notDE, trinotate, by.x="rn", by.y="#gene_id", all.x=T)
venom_500_notDE_sigp <- subset(venom_500_notDE_sigp, venom_500_notDE_sigp$SignalP=="Y")

#were they de vs other tissues
head <- fread("output/03_deseq2/tissue_itWT_LRT/Venom/venom_vs_tissue_res/Venom_vs_head_res.csv")
abdo <- fread("output/03_deseq2/tissue_itWT_LRT/Venom/venom_vs_tissue_res/Venom_vs_abdomen_res.csv")
thorax <- fread("output/03_deseq2/tissue_itWT_LRT/Venom/venom_vs_tissue_res/Venom_vs_thorax_res.csv")
ovary <- fread("output/03_deseq2/tissue_itWT_LRT/Venom/venom_vs_tissue_res/Venom_vs_ovaries_res.csv")
ovary_sig <- subset(ovary, padj<0.05)
# add columns for sig or not
venom_500_notDE_sigp$head_sig <- ifelse(venom_500_notDE_sigp$rn %in% head$rn, "TRUE", "FALSE")
venom_500_notDE_sigp$thorax_sig <- ifelse(venom_500_notDE_sigp$rn %in% thorax$rn, "TRUE", "FALSE")
venom_500_notDE_sigp$abdomen_sig <- ifelse(venom_500_notDE_sigp$rn %in% abdo$rn, "TRUE", "FALSE")
venom_500_notDE_sigp$ovary_sig <- ifelse(venom_500_notDE_sigp$rn %in% ovary_sig$rn, "TRUE", "FALSE")


# DE against all but ovary
de_abdo <- subset(venom_500_notDE_sigp, abdomen_sig==T)
de_head <- subset(venom_500_notDE_sigp, head_sig==T)
de_thorax <- subset(venom_500_notDE_sigp, thorax_sig==T)
de_all_others <- full_join(de_abdo, full_join(de_head, de_thorax))
de_all_others_res <- merge(de_all_others, ovary)
fwrite(de_all_others_res, "output/03_deseq2/tissue_itWT_LRT/Venom/venom_vs_tissue_res/de_all_but_ovary.csv")
