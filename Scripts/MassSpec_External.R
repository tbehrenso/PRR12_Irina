library(tidyverse)
library(ggplot2)
library(readxl)



# Example: read multiple replicates
GFP_PRR12 <- read_xlsx('/Users/tbehr/Desktop/MassSpecDataset_Intepreted.xlsx', sheet = 'GFP_PRR12')
CTRL1 <-  read_xlsx('/Users/tbehr/Desktop/MassSpecDataset_Intepreted.xlsx', sheet = 'Ctrl1')


prr12_df <- GFP_PRR12 %>%
  select(Accession, value = `# PSMs`) %>%
  group_by(Accession) %>%
  summarise(PRR12 = sum(value), .groups = "drop")

ctrl_df <- CTRL1 %>%
  select(Accession, value = `# PSMs`) %>%
  group_by(Accession) %>%
  summarise(CTRL = sum(value), .groups = "drop")

merged <- full_join(prr12_df, ctrl_df, by = "Accession")

merged <- merged %>%
  mutate(
    PRR12 = replace_na(PRR12, 0),
    CTRL  = replace_na(CTRL, 0)
  )

pseudo <- 1  # or 0.5

results <- merged %>%
  mutate(
    PRR12_adj = PRR12 + pseudo,
    CTRL_adj  = CTRL + pseudo,
    log2FC = log2(PRR12_adj / CTRL_adj)
  )


results <- results %>%
  arrange(desc(log2FC))

results <- results %>%
  mutate(
    mean_abundance = (PRR12_adj + CTRL_adj) / 2,
    log10_mean = log10(mean_abundance)
  )

ggplot(results, aes(x = log10_mean, y = log2FC)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(
    y = "log2 Fold Change (PRR12 / CTRL)",
    x = "log10 Mean Abundance"
  )


# Add gene symbol
library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

mapping <- getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters = "uniprotswissprot",
  values = results$Accession,
  mart = mart
)

mapping <- mapping %>%
  filter(hgnc_symbol != "") %>%
  distinct(uniprotswissprot, .keep_all = TRUE)

results <- results %>%
  left_join(mapping, by = c("Accession" = "uniprotswissprot")) %>%
  rename(GeneSymbol = hgnc_symbol)

write.csv(results, '/Users/tbehr/Desktop/MassSpecFoldChange.csv',quote = F, row.names = F)

