
library(GenomicRanges)
library(DESeq2)
library(ggplot2)

consensus_WT <- read.table('/Users/tbehr/Desktop/05_consensus_peaks/WT.seacr.consensus.peak_counts.bed', sep = '\t', header = F)
consensus_HET <- read.table('/Users/tbehr/Desktop/05_consensus_peaks/HET.seacr.consensus.peak_counts.bed', sep = '\t', header = F)
consensus_KO <- read.table('/Users/tbehr/Desktop/05_consensus_peaks/KO.seacr.consensus.peak_counts.bed', sep = '\t', header = F)

consensus_colnames <- c('chr','start','end','indiv_starts','indiv_ends',
                        'total_signal','max_signal','coords_max_signal','filenames','file_count')

colnames(consensus_WT) <- consensus_colnames
colnames(consensus_HET) <- consensus_colnames
colnames(consensus_KO) <- consensus_colnames




# load individual sample beds and evaluate overall peak strength
bed_colnames <- c('chr','start','end','total_signal','max_signal','coords_max_signal')
bed_WT_1 <- read.table('/Users/tbehr/Desktop/seacr/WT_R1.seacr.peaks.stringent.bed', sep='\t', header=F, col.names=bed_colnames)
bed_WT_2 <- read.table('/Users/tbehr/Desktop/seacr/WT_R2.seacr.peaks.stringent.bed', sep='\t', header=F, col.names=bed_colnames)
bed_WT_3 <- read.table('/Users/tbehr/Desktop/seacr/WT_R3.seacr.peaks.stringent.bed', sep='\t', header=F, col.names=bed_colnames)
bed_HET_1 <- read.table('/Users/tbehr/Desktop/seacr/HET_R1.seacr.peaks.stringent.bed', sep='\t', header=F, col.names=bed_colnames)
bed_HET_2 <- read.table('/Users/tbehr/Desktop/seacr/HET_R2.seacr.peaks.stringent.bed', sep='\t', header=F, col.names=bed_colnames)
bed_HET_3 <- read.table('/Users/tbehr/Desktop/seacr/HET_R3.seacr.peaks.stringent.bed', sep='\t', header=F, col.names=bed_colnames)
bed_KO_1 <- read.table('/Users/tbehr/Desktop/seacr/KO_R1.seacr.peaks.stringent.bed', sep='\t', header=F, col.names=bed_colnames)
bed_KO_2 <- read.table('/Users/tbehr/Desktop/seacr/KO_R2.seacr.peaks.stringent.bed', sep='\t', header=F, col.names=bed_colnames)
bed_KO_3 <- read.table('/Users/tbehr/Desktop/seacr/KO_R3.seacr.peaks.stringent.bed', sep='\t', header=F, col.names=bed_colnames)


bed_combined <- data.frame(
  value = c(bed_WT_2$max_signal,bed_WT_3$max_signal,
            bed_HET_1$max_signal,bed_HET_2$max_signal,bed_HET_3$max_signal,
            bed_KO_1$max_signal,bed_KO_3$max_signal),
  condition = rep(c('WT','HET','KO'),
                  times = c(length(bed_WT_2$max_signal) + length(bed_WT_3$max_signal),
                            length(bed_HET_1$max_signal) + length(bed_HET_2$max_signal) + length(bed_HET_3$max_signal),
                            length(bed_KO_1$max_signal) + length(bed_KO_3$max_signal))),
  replicate = factor(rep(rep(1:3, each = 1),
                         times = c(length(bed_WT_2$max_signal) + length(bed_WT_3$max_signal),
                                   length(bed_HET_1$max_signal) + length(bed_HET_2$max_signal) + length(bed_HET_3$max_signal),
                                   length(bed_KO_1$max_signal) + length(bed_KO_3$max_signal))))
)

ggplot(bed_combined, aes(x = value, fill = condition)) +
  geom_histogram(
    aes(y = after_stat(density)),
    alpha = 0.4,
    position = "identity",
    bins = 30
  ) +
  theme_minimal() +
  labs(
    title = "Overlayed Histogram of Value Distributions",
    x = "Value",
    y = "Count"
  ) + xlim(c(0, 50))












