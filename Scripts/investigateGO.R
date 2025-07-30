library(org.Hs.eg.db)
library(AnnotationDbi)

# Load Data
NEU_vst <- read.csv('data/NEU_heatmap_vst.csv', row.names = 1)
NPC_vst <- read.csv('data/NPC_heatmap_vst.csv', row.names = 1)

NEU_norm <- read.csv('data/NEU_counts_normalized.csv', row.names = 1)
NPC_norm <- read.csv('data/NPC_counts_normalized.csv', row.names = 1)




go_name <- 'signal transduction in response to DNA damage'
go_id <- "GO:0042770"
CELL_TYPE <- 'NEU'

entrez_ids <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = go_id,
                                    keytype = "GO",
                                    columns = c("ENTREZID", "SYMBOL"))

genes_of_go <- unique(entrez_ids$SYMBOL)
genes_for_heatmap <- genes_of_go

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
                         main = paste0(go_id, ': ', go_name),
                         cluster_rows=T, cluster_cols=FALSE,
                         column_names_side=c('top'), angle_col = c('45'),
                         row_names_side=c('left'),
                         heatmap_legend_param = list(title = NULL)
)
