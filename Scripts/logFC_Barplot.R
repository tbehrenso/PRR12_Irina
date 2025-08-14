library(ggplot2)

# Load Data
neu_ko_wt <- read.csv('Results/DEGs/DEGs_NEU_KO_vs_WT.csv', row.names = 1)

# Get Gene List
gene_lists <- read_xlsx('data/GeneLists.xlsx', sheet = 'Other')

# Choose gene subset
genes_for_heatmap <- unlist(as.data.frame(gene_lists[,c("CellCycle")]),use.names=F)
genes_for_heatmap <- genes_for_heatmap[!is.na(genes_for_heatmap)]

genes_for_heatmap <- c('CDK1','CDK2','CDK3','CDK6','CCNA2','CCNB1','CCNB2','CCND2','CCND2-AS1',
                       'CCNE2','CDKN1A','PKMYT1','CDC25A','CDC25C')

# Choose Cell Type
CELL_TYPE <- 'NEU'

gene_logfc_df <- data.frame(Gene = genes_for_heatmap, logfc = neu_ko_wt[genes_for_heatmap,]$log2FoldChange)
gene_logfc_df <- gene_logfc_df %>% mutate(regulation = ifelse(logfc > 0, 'up', 'down'))

ggplot(data=gene_logfc_df, aes(x=Gene, y=logfc, fill=regulation)) +
  geom_bar(stat="identity")


# --------------------------
# Alternative

GO_DEGs <- arrange(GO_DEGs, desc(log2FoldChange))

gene_logfc_df <- data.frame(Gene = rownames(GO_DEGs), logfc = GO_DEGs$log2FoldChange, padj = GO_DEGs$padj) %>% mutate(regulation = ifelse(logfc > 0, 'up', 'down'))

gene_logfc_df$Gene <- factor(gene_logfc_df$Gene, levels = gene_logfc_df$Gene)

gene_logfc_df$padj[gene_logfc_df$padj > 0.05] <- NA

ggplot(data=gene_logfc_df, aes(x=Gene, y=logfc, color=padj)) +
  geom_bar(stat="identity") +
  scale_colour_gradient(low='#bb39db',high='#3d2e47')
