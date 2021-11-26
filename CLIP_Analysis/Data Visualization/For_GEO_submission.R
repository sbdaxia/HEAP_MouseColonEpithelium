library(tidyverse)

peak_info <- readRDS("Datafiles/merged-peaks-mirs-200-12232019-withID.rds")
write_csv(as.data.frame(peak_info), "Peak_information.csv")

load("../Merged_Analysis/merged_peak_analysis.rda")
write.csv(as.data.frame(hfk.hf.res), "Peak_differential_analysis.csv", row.names = TRUE)
write.csv(as.data.frame(DESeq2::counts(dds.peaks.hf.hfk, normalize = TRUE)), "Peak_normalized_counts.csv", row.names = TRUE)


