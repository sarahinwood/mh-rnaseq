library(data.table)
library(dplyr)

head <- fread("output/deseq2/tissue_LRT/head/head_sp_LRT_annots.csv")
thorax <- fread("output/deseq2/tissue_LRT/thorax/thorax_sp_LRT_annots.csv")
abdo <- fread("output/deseq2/tissue_LRT/abdo/abdo_sp_LRT_annots.csv")
ovary <- fread("output/deseq2/tissue_LRT/ovary/ovary_sp_LRT_annots.csv")
venom <- fread("output/deseq2/tissue_LRT/venom/venom_sp_LRT_annots.csv")

##join all tables
all_tissue_degs <- full_join(head, thorax)
all_tissue_degs <- full_join(abdo, all_tissue_degs)
all_tissue_degs <- full_join(ovary, all_tissue_degs)
all_tissue_degs <- full_join(venom, all_tissue_degs)

##filter out those with no Blast annots
unann_degs <- subset(all_tissue_degs, sprot_Top_BLASTX_hit=="")
fwrite(list(unann_degs$rn), "output/deseq2/tissue_LRT/unann_degs_ids.txt")
