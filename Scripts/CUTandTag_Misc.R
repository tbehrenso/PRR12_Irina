


# ---------------------------------------------------------------
# Piechart showing which feature types (promoter, UTR, etc) the peaks fall on
# ---------------------------------------------------------------
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ChIPseeker)
library(ggplot2)
library(dplyr)

peaks_WT <- readPeakFile('/Users/tbehr/Desktop/KO.macs2.consensus.peak_counts.bed')
peaks_WT <- readPeakFile('/Users/tbehr/Desktop/peaks_macs2_broad_joined0_support2.bed')


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peak_anno <- annotatePeak(
  peaks_WT,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

feature_df <- peak_anno@annoStat

df_grouped <- feature_df %>%
  mutate(
    Feature_group = case_when(
      grepl("^Promoter", Feature) ~ "Promoter",
      grepl("Exon", Feature) ~ "Exon",
      grepl("Intron", Feature) ~ "Intron",
      Feature %in% c("5' UTR", "3' UTR") ~ Feature,
      Feature == "Downstream (<=300)" ~ "Downstream",
      Feature == "Distal Intergenic" ~ "Intergenic",
      TRUE ~ Feature
    )
  )
df_summary <- df_grouped %>%
  group_by(Feature_group) %>%
  summarise(Frequency = sum(Frequency)) %>%
  ungroup()



ggplot(df_summary, aes(x = "", y = Frequency, fill = Feature_group)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  theme_void() +
  labs(title = "Genomic Annotation of Peaks")



# ---------------------------------------------------------------
# Extract diffbind consensus peaks into a bed file
# ---------------------------------------------------------------
consensus_peaks <- dba.peakset(db, bRetrieve = TRUE)
rtracklayer::export(consensus_peaks, "/Users/tbehr/Desktop/diffbind_consensus_peaks.bed")


# ---------------------------------------------------------------
# From ChIP Atlas Enrichment of PRR12 Peaks, see over-represented TFs
# ---------------------------------------------------------------
library(dplyr)
library(tidyverse)

Prr12ChipAtlas_TFenrichment <- read.table('/Users/tbehr/Desktop/SanRaffaele/Projects/PRR12_Irina/data/External/Prr12ChipAtlas_TFenrichment.tsv', sep='\t')
colnames(Prr12ChipAtlas_TFenrichment) <- c('Dataset','name','Symbol','Source','CellType','PeakNum','PeakRatio','BackgroundRatio','pValue','qValue','FoldChange')

summary_df <- Prr12ChipAtlas_TFenrichment %>%
  mutate(neg_qValue = -qValue) %>% 
  filter(neg_qValue > 4, FoldChange > 5) %>%
  filter(!(Symbol %in% c('GFP','Epitope tags'))) %>% 
  group_by(Symbol) %>%
  summarise(
    n_hits = n(),
    mean_fold = mean(FoldChange, na.rm = TRUE),
    max_fold  = max(FoldChange, na.rm = TRUE),
    mean_qval = mean(neg_qValue, na.rm = TRUE),
    max_qval  = max(neg_qValue, na.rm = TRUE),
    .groups = "drop"
  )

summary_df <- summary_df %>%
  mutate(
    z_mean_fold = scale(log10(mean_fold)),
    z_max_fold  = scale(log10(max_fold)),
    z_mean_qval = scale(-mean_qval),
    z_max_qval  = scale(-max_qval),
    z_n_hits    = scale(n_hits)
  )

summary_df <- summary_df %>%
  mutate(
    composite_score =
      z_mean_fold +
#      z_max_fold +
      2*z_mean_qval +
#      z_max_qval +
      z_n_hits
  ) %>% 
  select(-one_of(c('z_mean_fold','z_max_fold','z_mean_qval','z_max_qval','z_n_hits'))) %>% 
  arrange(desc(composite_score))

# Alternate that uses rank scores instead of z-values
summary_df_alt <- Prr12ChipAtlas_TFenrichment %>%
  separate(PeakRatio, c('PeakOverlap','Prr12Peaks'),sep='/', convert = T) %>% 
  filter(CellType=='293') %>% 
  group_by(Symbol) %>%
  summarise(
    n_hits = n(),
    mean_fold = mean(FoldChange, na.rm = TRUE),
    max_fold  = max(FoldChange, na.rm = TRUE),
    mean_qval = mean(qValue, na.rm = TRUE),
    min_qval  = min(qValue, na.rm = TRUE),
    mean_percentOverlap = mean(PeakOverlap / Prr12Peaks),
    .groups = "drop"
  ) %>%
  mutate(
    rank_n_hits    = percent_rank(n_hits),
    rank_mean_fold = percent_rank(mean_fold),
    rank_max_fold  = percent_rank(max_fold),
    
    # lower q-value = better
    rank_mean_qval = percent_rank(-mean_qval),
    rank_min_qval  = percent_rank(-min_qval)
  ) %>%
  mutate(
    composite_score =
      (
#        rank_n_hits +
          rank_mean_fold +
          rank_mean_qval
      ) / 3
  ) %>%
  select(-one_of(c('rank_mean_fold','rank_max_fold','rank_mean_qval','rank_min_qval','rank_n_hits'))) %>% 
  arrange(desc(composite_score))


write.csv(summary_df_alt,'/Users/tbehr/Desktop/Prr12ChipAtlas_TF_customized.csv', row.names = F)





