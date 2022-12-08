library(data.table)
library(ggplot2)
library(viridis)
library(tidyverse)
library(scales)

# tpm data
salmon_tpm_file = 'output/03_deseq2/salmon_TPM.csv'
salmon_tpms <- fread(salmon_tpm_file)

# mean tpm for each tissue
salmon_mean_tpms <- data.table()
salmon_mean_tpms$gene <-salmon_tpms$rn

salmon_mean_tpms$Head <- rowMeans(subset(salmon_tpms, select = c(5,6,7,20)), na.rm = TRUE)
salmon_mean_tpms$Thorax <- rowMeans(subset(salmon_tpms, select = c(14,15,16)), na.rm = TRUE)
salmon_mean_tpms$Abdomen <- rowMeans(subset(salmon_tpms, select = c(2,3,4,19)), na.rm = TRUE)
salmon_mean_tpms$Venom <- rowMeans(subset(salmon_tpms, select = c(17,18)), na.rm = TRUE)
salmon_mean_tpms$Ovaries <- rowMeans(subset(salmon_tpms, select = c(8,9,10,21,22)), na.rm = TRUE)
salmon_mean_tpms$Pupa <- rowMeans(subset(salmon_tpms, select = c(11,12,13)), na.rm = TRUE)
fwrite(salmon_mean_tpms, "output/03_deseq2/salmon_TPMs_tissue_means.csv")

melted_mean_tpms <- melt(salmon_mean_tpms, value.name="TPM", variable.name="Tissue", id=c(1))
melted_mean_tpms_200 <- subset(melted_mean_tpms, TPM>200)

# plot
png("tpm.png", units="in", width=5, height=5, res=300)
ggplot(melted_mean_tpms, aes(y=TPM, x=Tissue))+
  geom_jitter(alpha=0.7)+
  scale_y_continuous(breaks = seq(0, 110000, by = 10000), labels = comma)+
  theme_bw()
dev.off()

# without larva or head max
ggplot(melted_mean_tpms, aes(y=TPM, x=Tissue))+
  geom_jitter(alpha=0.7)+
  theme_bw()+
  ylim(-10, 50000)

melted_tpms_above_10 <- subset(melted_mean_tpms, TPM>25 & TPM<1000)

#facetted by spread
ggplot(melted_mean_tpms, aes(x=TPM))+
  geom_density()+
  theme_bw()+
  scale_x_log10(breaks=c(0.01, 0.1, 1, 10, 100, 1000, 10000, 100000))+
  facet_wrap(~Tissue)

# coloured by tissue

melted_mean_tpms$Tissue <- factor(melted_mean_tpms$Tissue, levels = c("Pupa", "Head", "Thorax", "Abdomen", "Venom", "Ovaries"))
ggplot(melted_mean_tpms, aes(x=TPM, fill=Tissue))+
  geom_density(alpha=0.6)+
  theme_bw()+
  scale_x_log10(breaks=c(0.01, 0.1, 1, 10, 100, 1000, 10000, 100000))+
  scale_fill_viridis(discrete=T)+
  ylab("Density")
# labels on x axis are whole numbers already adjusted for log10 scaling I believe
  # most genes have less than 1 transcript per million



ggplot(melted_mean_tpms, aes(x=TPM))+
  geom_density(alpha=0.7)+
  theme_bw()+
  scale_x_continuous(trans='log10')+
  facet_wrap(~Tissue)
  
  #
  
melted_tpms <- melt(salmon_tpms, value.name="TPM", variable.name="Sample", id=c(1))
melted_tpms$Sample <- gsub("Mhyp_old_", "Mhold_", melted_tpms$Sample)
melted_tpms$tissue <- tstrsplit(melted_tpms$Sample, "_", keep=2)
melted_tpms$tissue <- tstrsplit(melted_tpms$tissue, "[[:digit:]]", keep=1)

# summary stats - mean TPM doesn't mean anything - they're adjusted to all add to 1mil, and each sample has same number of genes
group_by(melted_tpms, tissue) %>%
  summarise(
    count = n(),
    above5 = sum(TPM>5),
    above10 = sum(TPM>10),
    above100 = sum(TPM>100),
    above1000 = sum(TPM>1000),
    above10000 = sum(TPM>10000),
    sd = sd(TPM, na.rm = TRUE),
    median = median(TPM, na.rm = TRUE),
    IQR = IQR(TPM, na.rm = TRUE),
    max = max(TPM))

group_by(melted_mean_tpms, Tissue) %>%
  summarise(
    count = n(),
    above5 = sum(TPM>5),
    above10 = sum(TPM>10),
    above100 = sum(TPM>100),
    above200 = sum(TPM>200),
    above1000 = sum(TPM>1000),
    above10000 = sum(TPM>10000),
    sd = sd(TPM, na.rm = TRUE),
    median = median(TPM, na.rm = TRUE),
    IQR = IQR(TPM, na.rm = TRUE),
    max = max(TPM))

summary_chars <- capture.output(print(summary))

fwrite(summary_chars, "output/03_deseq2/TPM_stats.txt")



## example colouring points that are specifically upregulated in each tissue
# tissue specific DEGs
#Head <- fread('output/03_deseq2/tissue_itWT/Head/Head_sp_annots.csv')
#Head_up <- subset(Head, log2FoldChange > 0)
#head_mean_tpms <- salmon_mean_tpms
#head_mean_tpms$Tissue <- paste("Head")
#head_mean_tpms$TPM <- rowMeans(subset(salmon_tpms, select = c(5,6,7,20)), na.rm = TRUE)
#head_mean_tpms$tissue_specific <- ifelse(head_mean_tpms$gene %in% Head_up$rn, "Head-specific", "Non-specific")
#head_mean_tpms <- data.table(head_mean_tpms)
#ggplot(head_mean_tpms %>%
#         arrange(desc(tissue_specific)), aes(y=TPM, x=Tissue, color=tissue_specific))+
#  geom_jitter(alpha=0.7)+
#  theme_bw()+
#  scale_colour_viridis(discrete=T, direction=-1)+
#  ylim(-10, 105000)
