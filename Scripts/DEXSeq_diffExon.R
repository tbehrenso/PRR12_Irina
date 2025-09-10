library(DEXSeq)
library(tidyverse)
library(ggplot2)

PRR12_ID <- 'ENSG00000126464'
CELL_TYPE <- 'NPC'

# ------------------------------------------
# Preparation

# Prepare Annotation

# Count Reads

# ------------------------------------------
# Read in Prepared Annotation and Counts

sample_info <- read.csv('data/sample_info.csv')
sample_info$CellType <- factor(sample_info$cellType)
sample_info$condition <- factor(sample_info$condition)

count_files <- c('data/HTSeq_Counts/NEU_WT_1.count.txt','data/HTSeq_Counts/NEU_WT_2.count.txt','data/HTSeq_Counts/NEU_WT_3.count.txt','data/HTSeq_Counts/NEU_WT_4.count.txt',
                 'data/HTSeq_Counts/NEU_HET_1.count.txt','data/HTSeq_Counts/NEU_HET_2.count.txt','data/HTSeq_Counts/NEU_HET_3.count.txt','data/HTSeq_Counts/NEU_HET_4.count.txt',
                 'data/HTSeq_Counts/NEU_KO_1.count.txt','data/HTSeq_Counts/NEU_KO_2.count.txt','data/HTSeq_Counts/NEU_KO_3.count.txt','data/HTSeq_Counts/NEU_KO_4.count.txt',
                 'data/HTSeq_Counts/NPC_WT_1.count.txt','data/HTSeq_Counts/NPC_WT_2.count.txt','data/HTSeq_Counts/NPC_WT_3.count.txt','data/HTSeq_Counts/NPC_WT_4.count.txt',
                 'data/HTSeq_Counts/NPC_HET_1.count.txt','data/HTSeq_Counts/NPC_HET_2.count.txt','data/HTSeq_Counts/NPC_HET_3.count.txt','data/HTSeq_Counts/NPC_HET_4.count.txt',
                 'data/HTSeq_Counts/NPC_KO_1.count.txt','data/HTSeq_Counts/NPC_KO_2.count.txt','data/HTSeq_Counts/NPC_KO_3.count.txt','data/HTSeq_Counts/NPC_KO_4.count.txt')

if(CELL_TYPE == 'NEU'){
  sample_info <- sample_info[1:12,]
  count_files <- count_files[1:12]
} else if(CELL_TYPE == 'NPC'){
  sample_info <- sample_info[13:24,]
  count_files <- count_files[13:24]
}


# ------------------------------------------
# Standard Analysis Workflow

# Build DEXSeqDataSet
dxd <- DEXSeqDataSetFromHTSeq(count_files, 
                              sampleData = sample_info,
                              design = ~ sampleName + exon + condition:exon,
                              flattenedfile = 'data/ensembl.h37.87.DEXSeq.chr.gff')

formulaFullModel    =  ~ sampleName + exon + condition:exon
formulaReducedModel =  ~ sampleName + exon

# Processing
dxd = estimateSizeFactors(dxd)

gene_counts <- rowsum(featureCounts(dxd), group = geneIDs(dxd))
keep_genes <- rowSums(gene_counts >= 2000) >= 8

genes_to_keep <- names(keep_genes[keep_genes])
keep <- geneIDs(dxd) %in% genes_to_keep

dxd_filtered <- dxd[keep,]

# dxd_filtered <- dxd[geneIDs(dxd)==PRR12_ID,]

dxd_filtered = estimateDispersions(dxd_filtered)

###### RData file saved from HERE (after estimateDispersions, before diff. exon usage calc)

plotDispEsts(dxd_filtered)

# Test for Differential Exon Usage
dxd_filtered <- testForDEU(dxd_filtered)
dxd_filtered <- estimateExonFoldChanges(dxd_filtered, fitExpToVar = "condition")
results <- DEXSeqResults(dxd_filtered)

###### _diffExp.RData file saved from HERE (after estimateDispersions, before diff. exon usage calc)

plotDEXSeq(results, PRR12_ID, legend=TRUE, displayTranscripts = T, cex.axis=1.2, cex=1.3, lwd=2)
plotDEXSeq(results, PRR12_ID, legend=TRUE, expression=F, norCounts=T, displayTranscripts = T, cex.axis=1.2, cex=1.3, lwd=2)






# ------------------------------------------------------------------------------------
#  BELOW is just copied from Simone's DEXSeq Analysis (obviously needs modification before use here)
# ------------------------------------------------------------------------------------


dxd <- testForDEU(dxd)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition")
results <- DEXSeqResults(dxd)

# count differentially expressed EXONS
table(results$padj < 0.05 )
# count differentially affected genes
table(tapply(results$padj < 0.1, results$groupID, any))


# filtered_results <- results[!is.na(results$padj) & !is.infinite(results$log2fold_PW1_EDO), ]

# plot a specific gene
plotDEXSeq(results, "ENSG00000224078.15", legend=TRUE, displayTranscripts = T, cex.axis=1.2, cex=1.3, lwd=2)


results_df <- results %>% 
  as.data.frame %>%
  mutate(logp=-log10(padj)) %>% 
  mutate(highlight=logp > -log10(0.05) & (log2fold_PW1_CTRL < -1 | log2fold_PW1_CTRL > 1))

# Volcano Plot
ggplot(as.data.frame(results_df), aes(x=log2fold_PW1_CTRL, y=logp, color=highlight)) +
  geom_point(size=2) +
  theme_bw() +
  ylab('-log10(adjusted p-value)') + xlab('log2(fold-change)') +
  theme(text=element_text(size=20))
#scale_color_manual(values=c('#7f7f7f','#cd1e05'))
