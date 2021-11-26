library(data.table)
library(tidyverse)
library(rtracklayer)
library(Rsamtools)
library(DESeq2)
library(Rsubread)
library(CLIPanalyze)
library(RColorBrewer)
library(gplots)
library(Biostrings)

peak.data.hf<- readRDS("peakdata.HF.rds")
peaks.all.hf <- peak.data.hf$peaks
peaks.all.hf <- subset(peaks.all.hf, width > 20)
peaks.all.hf <- subset(peaks.all.hf, log2FC > 0)
count.threshold <- 10
norm.counts <- counts(peak.data.hf$peak.counts, normalized = TRUE)
norm.counts <- norm.counts[names(peaks.all.hf), ]
colnames(norm.counts)[1:6] <- c(paste0("HF", 1:3), paste0("HF-IC", 1:3))
selected.peaks <- rowMeans(norm.counts[, paste0("HF", 1:3)]) > count.threshold 
peaks.all.hf <- peaks.all.hf[selected.peaks, ]
padj.threshold <- 5*1e-2
hf.peaks <- subset(peaks.all.hf, padj < padj.threshold)
length(hf.peaks)
table(hf.peaks$annot)

peak.data.hfk<- readRDS("peakdata.HFK.rds")
peaks.all.hfk <- peak.data.hfk$peaks
peaks.all.hfk <- subset(peaks.all.hfk, width > 20)
peaks.all.hfk <- subset(peaks.all.hfk, log2FC > 0)
norm.counts <- counts(peak.data.hfk$peak.counts, normalized = TRUE)
norm.counts <- norm.counts[names(peaks.all.hfk), ]
colnames(norm.counts)[1:6] <- c(paste0("HFK", 1:3), paste0("HFK-IC", 1:3))
selected.peaks <- rowMeans(norm.counts[, paste0("HFK", 1:3)]) > count.threshold 
peaks.all.hfk <- peaks.all.hfk[selected.peaks, ]
hfk.peaks <- subset(peaks.all.hfk, padj < padj.threshold)
length(hfk.peaks)
table(hfk.peaks$annot)


# Make staked bargraphs for peak fractions
library(ggplot2)
library(tidyverse)
library(reshape2)
hf_peak_frac <- read_csv("feature_frac_hf_peaks.csv")
hfk_peak_frac <- read_csv("feature_frac_hfk_peaks.csv")
hf_peak_frac_df <- melt(hf_peak_frac)
colnames(hf_peak_frac_df) <- c("region","rank","fraction")
hf_peak_frac_df <- hf_peak_frac_df %>% filter(!region %in% c("utr3*","utr5*"))
hfk_peak_frac_df <- melt(hfk_peak_frac)
colnames(hfk_peak_frac_df) <- c("region","rank","fraction")
hfk_peak_frac_df <- hfk_peak_frac_df %>% filter(!region %in% c("utr3*","utr5*"))

pdf("PDF_figure/HF_peak_genomic_fraction.pdf",
    height = 4,
    width = 6)
ggplot(hf_peak_frac_df, aes(fill=region, y=fraction, x=rank)) + 
  geom_bar(position="fill", stat="identity")+
  xlab("Top n peaks ranked by adjusted p-vaue") +
  ylab("Distribution of peaks\nacross genomic annotations") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size = 0.5)) + 
  scale_y_continuous(expand = c(0, 0))
dev.off()

pdf("PDF_figure/HFK_peak_genomic_fraction.pdf",
    height = 4,
    width = 6)
ggplot(hfk_peak_frac_df, aes(fill=region, y=fraction, x=rank)) + 
  geom_bar(position="fill", stat="identity")+
  xlab("Top n peaks ranked by adjusted p-vaue") +
  ylab("Distribution of peaks\nacross genomic annotations") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size = 0.5)) + 
  scale_y_continuous(expand = c(0, 0))
dev.off()