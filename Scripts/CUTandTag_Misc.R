


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















