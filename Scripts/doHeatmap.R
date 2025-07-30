library(pheatmap)
library(tidyverse)
library(dplyr)
library(readxl)

# Load Data
NEU_vst <- read.csv('data/NEU_heatmap_vst.csv', row.names = 1)
NPC_vst <- read.csv('data/NPC_heatmap_vst.csv', row.names = 1)

NEU_norm <- read.csv('data/NEU_counts_normalized.csv', row.names = 1)
NPC_norm <- read.csv('data/NPC_counts_normalized.csv', row.names = 1)

# Get Gene List
gene_lists <- read_xlsx('data/GeneLists.xlsx', sheet = 'Other')

# Choose gene subset
genes_for_heatmap <- unlist(as.data.frame(gene_lists[,c("CellCycle")]),use.names=F)
genes_for_heatmap <- genes_for_heatmap[!is.na(genes_for_heatmap)]

# Search for specific gene
NEU_norm[grep('CIT',rownames(NEU_norm)),]

# Choose Cell Type
CELL_TYPE <- 'NEU'

# Subset Matrix
if(CELL_TYPE=='NPC'){
  genes_for_heatmap <- genes_for_heatmap[genes_for_heatmap %in% rownames(NPC_norm)]
  df_subset <- NPC_norm[genes_for_heatmap,]
} else if (CELL_TYPE=='NEU'){
  genes_for_heatmap <- genes_for_heatmap[genes_for_heatmap %in% rownames(NEU_norm)]
  df_subset <- NEU_norm[genes_for_heatmap,]
}


ComplexHeatmap::pheatmap(as.matrix(df_subset),
                         scale = "row",
                         color = rev(RColorBrewer::brewer.pal(n =4, name = "RdYlBu")),
                         #color = inferno(256),
                         #color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                         main = "CellCycle",
                         cluster_rows=T, cluster_cols=FALSE,
                         column_names_side=c('top'), angle_col = c('45'),
                         row_names_side=c('left'),
                         heatmap_legend_param = list(title = NULL)
)

# Heatmap with ggplot
dea_counts_subset_scaled <- t(scale(t(dea_counts_subset)))
dea_counts_subset_modified <- data.frame(Gene=rownames(dea_counts_subset_scaled),dea_counts_subset_scaled)
dea_counts_subset_long <- pivot_longer(dea_counts_subset_modified, cols = -Gene)
dea_counts_subset_long$Gene <- factor(dea_counts_subset_long$Gene, levels = rev(features_of_interest))  # set genes as factor to order the rows

dea_counts_subset_long$name <- factor(dea_counts_subset_long$name, levels=c('WT_1','WT_2','WT_3','SNCA_1','SNCA_2','SNCA_3'))

ggplot(dea_counts_subset_long, aes(x=name, y=Gene, fill=value)) +
  geom_tile(width=0.99, height=0.99) +
  #scale_fill_viridis(option='rocket', discrete=F) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme_minimal() +
  theme(text=element_text(size=30), legend.title=element_blank(), legend.key.height = unit(2, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_discrete(position = "top") +
  labs(x=NULL,y=NULL)




