library(csaw)
library(dplyr)

HISTONE_MODIFICATION <- 'H3k27me3'
HISTONE_MODIFICATION <- 'H3k27ac'
HISTONE_MODIFICATION <- 'H3k9me3'


if(HISTONE_MODIFICATION == 'H3k27me3'){
  merged <- readRDS('/Users/tbehr/Desktop/SanRaffaele/Projects/PRR12_Irina/Results/CUTandTag/H3k27me3/csaw/merged_WTvKO_v5.RData')
  clustered <- readRDS('/Users/tbehr/Desktop/SanRaffaele/Projects/PRR12_Irina/Results/CUTandTag/H3k27me3/csaw/clustered_WTvKO_v5.RData')
  mean_cpm <- readRDS('/Users/tbehr/Desktop/SanRaffaele/Projects/PRR12_Irina/Results/CUTandTag/H3k27me3/csaw/mean_cpm_obj_v5.RData')
  directory <- '/Users/tbehr/Desktop/SanRaffaele/Projects/PRR12_Irina/Results/CUTandTag/H3k27me3/csaw'
} else if(HISTONE_MODIFICATION == 'H3k27ac'){
  merged <- readRDS('/Users/tbehr/Desktop/SanRaffaele/Projects/PRR12_Irina/Results/CUTandTag/H3k27ac/csaw/merged_WTvKO_v5.RData')
  clustered <- readRDS('/Users/tbehr/Desktop/SanRaffaele/Projects/PRR12_Irina/Results/CUTandTag/H3k27ac/csaw/clustered_WTvKO_v5.RData')
  mean_cpm <- readRDS('/Users/tbehr/Desktop/SanRaffaele/Projects/PRR12_Irina/Results/CUTandTag/H3k27ac/csaw/mean_cpm_obj_v5.RData')
  directory <- '/Users/tbehr/Desktop/SanRaffaele/Projects/PRR12_Irina/Results/CUTandTag/H3k27ac/csaw'
} else if(HISTONE_MODIFICATION == 'H3k9me3'){
  merged <- readRDS('/Users/tbehr/Desktop/SanRaffaele/Projects/PRR12_Irina/Results/CUTandTag/H3k9me3/csaw/merged_WTvKO_v6.RData')
  clustered <- readRDS('/Users/tbehr/Desktop/SanRaffaele/Projects/PRR12_Irina/Results/CUTandTag/H3k9me3/csaw/clustered_WTvKO_v6.RData')
  mean_cpm <- readRDS('/Users/tbehr/Desktop/SanRaffaele/Projects/PRR12_Irina/Results/CUTandTag/H3k9me3/csaw/mean_cpm_obj_WTvKO_v6.RData')
  directory <- '/Users/tbehr/Desktop/SanRaffaele/Projects/PRR12_Irina/Results/CUTandTag/H3k9me3/csaw'
}

merged$regions


tabcom <- merged$combined
is.sig <- tabcom$FDR <= 0.05
summary(is.sig)

table(tabcom$direction[is.sig])


tabbest <- merged$best
is.sig.pos <- (tabbest$rep.logFC > 0)[is.sig]
summary(is.sig.pos)

library(rtracklayer)
merged_sig_regions <- merged$regions[is.sig]

cpm_filter <- case_when(
  HISTONE_MODIFICATION == 'H3k27me3' ~ 2,
  HISTONE_MODIFICATION == 'H3k27ac' ~ 0.5,
  HISTONE_MODIFICATION == 'H3k9me3' ~ 1,
)

clustered_cpmFitlered <- clustered$regions[mean_cpm > cpm_filter]
clustered_cpmFitlered_stats <- clustered$stats[mean_cpm > cpm_filter,]


if(HISTONE_MODIFICATION == 'H3k27me3'){
  export(merged_sig_regions, paste0(directory,'/csaw_sig_merged_k27me3_v5.bed'))
  export(clustered$regions, paste0(directory,'/csaw_sig_clustered_k27me3_v5.bed'))
  export(clustered_cpmFitlered, paste0(directory,'/csaw_sig_clustered_cpmFilt_k27me3_v5.bed'))
  write.csv(clustered_cpmFitlered_stats, paste0(directory,'/csaw_sig_clustered_cpmFilt_k27me3_STATS_v5.csv'),quote = F, col.names = T)
} else if(HISTONE_MODIFICATION == 'H3k27ac'){
  export(merged_sig_regions, paste0(directory,'/csaw_sig_merged_k27ac_v5.bed'))
  export(clustered$regions, paste0(directory,'/csaw_sig_clustered_k27ac_v5.bed'))
  export(clustered_cpmFitlered, paste0(directory,'/csaw_sig_clustered_cpmFilt_k27ac_v5.bed'))
  write.csv(clustered_cpmFitlered_stats, paste0(directory,'/csaw_sig_clustered_cpmFilt_k27ac_STATS_v5.csv'),quote = F, col.names = T)
} else if(HISTONE_MODIFICATION == 'H3k9me3'){
  export(merged_sig_regions, paste0(directory,'/csaw_sig_merged_k9me3_v6.bed'))
  export(clustered$regions, paste0(directory,'/csaw_sig_clustered_k9me3_v6.bed'))
  export(clustered_cpmFitlered, paste0(directory,'/csaw_sig_clustered_cpmFilt_k9me3_v6.bed'))
  write.csv(clustered_cpmFitlered_stats, paste0(directory,'/csaw_sig_clustered_cpmFilt_k9me3_STATS_v6.csv'),quote = F, col.names = T)
}
















