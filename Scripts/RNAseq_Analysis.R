library(DESeq2)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(pheatmap)

# Load Data
raw_count_import <- read.delim('data/irina.all.counts.trimmed.tsv')

# Create count matrix
raw_count_matrix <- raw_count_import |> 
  select(paste(rep(c('NEU','NPC'),each=12), rep(c('WT','HET','KO'),each=4), c(1:4),sep = '_')) |> 
  as.matrix()
rownames(raw_count_matrix) <- raw_count_import$Geneid
  
# Sample info
sample_info <- read.csv('data/sample_info.csv')
sample_info$CellType <- factor(sample_info$CellType)
sample_info$Condition <- factor(sample_info$Condition)

# Remove genes with very low counts (prefiltering)
#raw_count_matrix <- raw_count_matrix[rowSums(raw_count_matrix) >= 5,]
raw_count_matrix <- raw_count_matrix[rowSums(raw_count_matrix >= 10) >= 4,]

# -----------------------------------------
# Visualize Count Data
# -----------------------------------------

logcounts <- as.matrix(log2(raw_count_matrix + 1))

par(cex.axis=0.8, mar=c(6, 4, 2, 2))
boxplot(logcounts, las=2, ylab="Log2(Counts)")
abline(h=median(logcounts), col="blue")

plot(rowMeans(logcounts), rowSds(logcounts), 
     main='Log2 Counts: sd vs mean')

# PCA
library(ggfortify)
library(ggrepel)

rlogcounts <- rlog(raw_count_matrix) # rlog from DEseq2
pcDat <- prcomp(t(rlogcounts))
autoplot(pcDat, data = sample_info, colour = 'CellType', shape = 'Condition', size=3) +
  geom_text_repel(aes(x=PC1, y=PC2, label=SampleName), box.padding = 0.8)

# Mean vs Variance
counts <- raw_count_matrix[,21:24]

gene_means <- rowMeans(counts)
gene_vars <- apply(counts, 1, var)

mean_var_df <- data.frame(
  mean = gene_means,
  variance = gene_vars
)

ggplot(mean_var_df, aes(x = mean, y = variance)) +
  geom_point(alpha = 0.4) +
  scale_x_log10() + 
  scale_y_log10() + 
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Mean Expression (log10)",
    y = "Variance (log10)",
    title = "Mean vs Variance of Gene Expression"
  )

# -----------------------------------------
# DEG analysis
# -----------------------------------------

CELL_TYPE <- 'NEU'

# Create DESeq object
if(CELL_TYPE == 'BOTH'){
  dds <- DESeqDataSetFromMatrix(countData = raw_count_matrix,
                                colData = sample_info,
                                design = ~ Condition + CellType)
  dds$Condition <- factor(dds$Condition, levels = c('WT','HET','KO'))
  dds <- DESeq(dds)
  res_HET_WT <- results(dds, contrast=c('Condition','HET','WT'))
  res_KO_WT <- results(dds, contrast=c('Condition','KO','WT'))
  
} else {
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(select(as.data.frame(raw_count_matrix), contains(CELL_TYPE))),
                                colData = sample_info[`if`(CELL_TYPE=='NEU', c(1:12), c(13:24)),],
                                design = ~ Condition)
  dds$Condition <- factor(dds$Condition, levels = c('WT','HET','KO'))
  #dds$condition <- relevel(dds$condition, ref = "WT")
  
  dds <- DESeq(dds)
  res_HET_WT <- results(dds, contrast=c('Condition','HET','WT'))
  res_KO_WT <- results(dds, contrast=c('Condition','KO','WT'))
}


# Visualization
library(EnhancedVolcano)

plotMA(res_KO_WT, ylim=c(-2,2))

EnhancedVolcano(res_KO_WT,
                lab = rownames(res_KO_WT),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = paste0(CELL_TYPE, ': KO vs WT'),
                subtitle = NULL)

# Plot gene
GENE_OF_INTEREST <- 'CDC25C'

plot_data <- plotCounts(dds, gene=GENE_OF_INTEREST, intgroup=c("Condition",'CellType'), 
                returnData=TRUE)

ggplot(plot_data, aes(x=Condition, y=count, color=CellType)) + 
  geom_point(position=position_jitter(w=0.1,h=0), size=4) +
  ylim(c(0, max(plot_data$count) + max(plot_data$count)*0.1)) +
  ggtitle(paste0(GENE_OF_INTEREST, ' expression')) +
  theme(plot.title = element_text(size = 30))

#
vsd <- vst(dds, blind=FALSE)
library("RColorBrewer")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$SampleName
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


# ----------------------------------------
# GO term
# -----------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Load results and get gene symbols
res_NEU_HET_vs_WT <- read.csv('Results/DEGs/DEGs_NEU_HET_vs_WT.csv')
res_NEU_KO_vs_WT <- read.csv('Results/DEGs/DEGs_NEU_KO_vs_WT.csv')
res_NPC_HET_vs_WT <- read.csv('Results/DEGs/DEGs_NPC_HET_vs_WT.csv')
res_NPC_KO_vs_WT <- read.csv('Results/DEGs/DEGs_NPC_KO_vs_WT.csv')


CELL_TYPE <- 'NEU'
MUTANT_TYPE <- 'KO'
REG_DIRECTION <- 'UP'


res_of_interest <- case_when(
  CELL_TYPE == 'NEU' & MUTANT_TYPE == 'HET' ~ res_NEU_HET_vs_WT,
  CELL_TYPE == 'NEU' & MUTANT_TYPE == 'KO' ~ res_NEU_KO_vs_WT,
  CELL_TYPE == 'NPC' & MUTANT_TYPE == 'HET' ~ res_NPC_HET_vs_WT,
  CELL_TYPE == 'NPC' & MUTANT_TYPE == 'KO' ~ res_NPC_KO_vs_WT)


gene_symbols <- if(REG_DIRECTION == 'UP'){
  res_of_interest$X[res_of_interest$padj<0.05 & res_of_interest$log2FoldChange > 1]
} else if (REG_DIRECTION == 'DOWN'){
  res_of_interest$X[res_of_interest$padj<0.05 & res_of_interest$log2FoldChange < -1]
} else {
  print('Invalid Regulation Direction')
}

# Convert gene symbols to Entrez IDs
gene_df <- bitr(gene_symbols, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)


# Remove duplicates (if any)
gene_entrez <- unique(gene_df$ENTREZID)
length(gene_entrez)

# ------------------------------
# GO Enrichment Analysis (Biological Process)
go_bp <- enrichGO(gene = gene_entrez,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  keyType = "ENTREZID",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

dotplot(go_bp, showCategory = 20) + ggtitle(paste0(CELL_TYPE, ' - ', MUTANT_TYPE, ' vs WT: ', REG_DIRECTION, 'regulated GO terms'))

# Extract terms
GO_TERM_OF_INTEREST <- 'DNA repair'

go_term_index <- which(go_bp@result$Description == GO_TERM_OF_INTEREST)

degs_of_go_term <- unlist(strsplit(go_bp@result$geneID[go_term_index], split = "/"))

writeClipboard(degs_of_go_term)

# Search for terms
go_term_search_indeces <- grep('damage', go_bp@result$Description)
go_bp@result$Description[go_term_search_indeces]


# ------------------------------
# KEGG Pathway Enrichment
kegg_res <- enrichKEGG(gene = gene_entrez,
                       organism = 'hsa',
                       keyType = "kegg",
                       pvalueCutoff = 0.05)
# Convert Entrez back to gene symbols for readability
kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

barplot(kegg_res, showCategory = 20) + ggtitle(paste0(CELL_TYPE, ' - ', MUTANT_TYPE, ' vs WT: ', REG_DIRECTION, 'regulated KEGG Pathways'))


# --------------------------------
# GSEA 
# --------------------------------

logfc_list <- res_of_interest$log2FoldChange
names(logfc_list) <- res_of_interest$X
logfc_list <- sort(logfc_list, decreasing=T)

# Run GSEA
gse <- gseGO(geneList = logfc_list,
             ont = 'BP',
             OrgDb = 'org.Hs.eg.db',
             keyType = 'SYMBOL'
)

dotplot(gse, showCategory=15, split=".sign") +  
  theme(text=element_text(size=24), axis.text.x = element_text(size=24), axis.text.y = element_text(size=24), axis.title.x = element_text(size=24)) + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=80)) +
  scale_size_area(max_size = 10) + 
  ggtitle(paste0('GSEA ', CELL_TYPE, ': ', MUTANT_TYPE, ' vs WT ')) + facet_grid(.~.sign)


emapplot(gse, showCategory = 10)

cnetplot(gse, categorySize="pvalue", foldChange=logfc_list, showCategory = 3)

ridgeplot(gse) + labs(x = "enrichment distribution")

gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)




















