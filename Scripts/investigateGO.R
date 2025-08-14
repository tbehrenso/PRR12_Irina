library(org.Hs.eg.db)
library(AnnotationDbi)

# Load Data
NEU_vst <- read.csv('data/NEU_heatmap_vst.csv', row.names = 1)
NPC_vst <- read.csv('data/NPC_heatmap_vst.csv', row.names = 1)

NEU_norm <- read.csv('data/NEU_counts_normalized.csv', row.names = 1)
NPC_norm <- read.csv('data/NPC_counts_normalized.csv', row.names = 1)




go_name <- 'Double-strand break repair'
go_id <- "GO:0043066"
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


# -----------------------------------------------
# Extract Significant Genes
# -----------------------------------------------

# Load DEGs
DEGS_KO_vs_WT <- read.csv('Results/DEGs/DEGs_NEU_KO_vs_WT.csv', row.names = 1)

GO_DEGs <- DEGS_KO_vs_WT[genes_of_go,c('log2FoldChange','padj')]

sig_DEGs <- GO_DEGs[abs(GO_DEGs$log2FoldChange) > 0.9 & GO_DEGs$padj < 0.06,]
sig_DEGs <- na.exclude(sig_DEGs)


writeClipboard(as.character(rownames(sig_DEGs)))
writeClipboard(as.character(sig_DEGs$log2FoldChange))
writeClipboard(as.character(sig_DEGs$padj))




