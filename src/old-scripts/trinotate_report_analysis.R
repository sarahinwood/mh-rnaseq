library(data.table)
library(trinotateR)

trinotate_annotation_report <- read_trinotate("data/trinotate_annotation_report.txt")
trinotate_summary <- summary_trinotate(trinotate_annotation_report)
write.csv(x = trinotate_summary, file = "output/trinotateR_summary.csv")
write.csv(x = trinotate_annotation_report, file = "output/annotation_report.csv")
