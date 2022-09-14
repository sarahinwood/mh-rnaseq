library(RnaSeqSampleSize)
library(DESeq2)

mh_dds <- readRDS("output/deseq2/mh_dds.rds")
counts <- data.frame(counts(mh_dds, normalized=F))

## stage
stage_counts <- counts[,c(10:12, 1:9, 13:21)]
stage_disp <- est_count_dispersion(stage_counts, group=c(rep(0,3), rep(1,18)))
stage_power_1 <- est_power_distribution(n=21, f=0.05, rho=1, distributionObject=stage_disp, repNumber=10000)
stage_power_2 <- est_power_distribution(n=21, f=0.05, rho=2, distributionObject=stage_disp, repNumber=10000)
stage_power_5 <- est_power_distribution(n=21, f=0.05, rho=5, distributionObject=stage_disp, repNumber=10000)
stage_power_10 <- est_power_distribution(n=21, f=0.05, rho=10, distributionObject=stage_disp, repNumber=10000)

# head
head_counts <- counts[,c(4,5,6,19,1,2,3,7:18,20,21)]
head_disp <- est_count_dispersion(stage_counts, group=c(rep(0,4), rep(1,17)))
head_power_1 <- est_power_distribution(n=21, f=0.05, rho=1, distributionObject=head_disp, repNumber=10000)
head_power_2 <- est_power_distribution(n=21, f=0.05, rho=2, distributionObject=head_disp, repNumber=10000)
head_power_5 <- est_power_distribution(n=21, f=0.05, rho=5, distributionObject=head_disp, repNumber=10000)
head_power_10 <- est_power_distribution(n=21, f=0.05, rho=10, distributionObject=head_disp, repNumber=10000)

# thorax
thorax_counts <- counts[,c(13,14,15,1:12,16:21)]
thorax_disp <- est_count_dispersion(stage_counts, group=c(rep(0,3), rep(1,18)))
thorax_power_1 <- est_power_distribution(n=21, f=0.05, rho=1, distributionObject=thorax_disp, repNumber=10000)
thorax_power_2 <- est_power_distribution(n=21, f=0.05, rho=2, distributionObject=thorax_disp, repNumber=10000)
thorax_power_5 <- est_power_distribution(n=21, f=0.05, rho=5, distributionObject=thorax_disp, repNumber=10000)
thorax_power_10 <- est_power_distribution(n=21, f=0.05, rho=10, distributionObject=thorax_disp, repNumber=10000)

# abdo
abdo_counts <- counts[,c(1,2,3,18,4:17,19:21)]
abdo_disp <- est_count_dispersion(abdo_counts, group=c(rep(0,4), rep(1,17)))
abdo_power_1 <- est_power_distribution(n=21, f=0.05, rho=1, distributionObject=abdo_disp, repNumber=10000)
abdo_power_2 <- est_power_distribution(n=21, f=0.05, rho=2, distributionObject=abdo_disp, repNumber=10000)
abdo_power_5 <- est_power_distribution(n=21, f=0.05, rho=5, distributionObject=abdo_disp, repNumber=10000)
abdo_power_10 <- est_power_distribution(n=21, f=0.05, rho=10, distributionObject=abdo_disp, repNumber=10000)

# ovaries
ov_counts <- counts[,c(7,8,9,20,21,1:6,10:19)]
ov_disp <- est_count_dispersion(stage_counts, group=c(rep(0,5), rep(1,16)))
ov_power_1 <- est_power_distribution(n=21, f=0.05, rho=1, distributionObject=ov_disp, repNumber=10000)
ov_power_2 <- est_power_distribution(n=21, f=0.05, rho=2, distributionObject=ov_disp, repNumber=10000)
ov_power_5 <- est_power_distribution(n=21, f=0.05, rho=5, distributionObject=ov_disp, repNumber=10000)
ov_power_10 <- est_power_distribution(n=21, f=0.05, rho=10, distributionObject=ov_disp, repNumber=10000)

# venom
ven_counts <- counts[,c(16,17,1:15,18:21)]
ven_disp <- est_count_dispersion(stage_counts, group=c(rep(0,2), rep(1,19)))
ven_power_1 <- est_power_distribution(n=21, f=0.05, rho=1, distributionObject=ov_disp, repNumber=10000)
ven_power_2 <- est_power_distribution(n=21, f=0.05, rho=2, distributionObject=ov_disp, repNumber=10000)
ven_power_5 <- est_power_distribution(n=21, f=0.05, rho=5, distributionObject=ov_disp, repNumber=10000)
ven_power_10 <- est_power_distribution(n=21, f=0.05, rho=10, distributionObject=ov_disp, repNumber=10000)


