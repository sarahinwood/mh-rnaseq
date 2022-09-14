library(data.table)
library(dplyr)
library(ggplot2)

venom_phyla <- fread("output/deseq2/tissue_itWT_LRT/venom/bacterial_phyla.csv")
venom_phyla_counts <- count(venom_phyla, phyla)
venom_phyla_counts <- subset(venom_phyla_counts, !(phyla=="NA"))

ggplot(venom_phyla_counts, aes(x=reorder(phyla, -n), y=n))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 65, hjust = 1, face = "italic"))+
  geom_col(alpha=0.8, fill="#440154FF", colour="#440154FF", width=0.9)+
  xlab("Bacterial Phyla")+
  ylab("Number of DEGs")
