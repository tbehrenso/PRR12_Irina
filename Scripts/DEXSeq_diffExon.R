library(DEXSeq)
library(tidyverse)
library(ggplot2)


# ------------------------------------------
# Preparation

# Prepare Annotation

# Count Reads

# Build DEXSeqDataSet
dxd <- DEXSeqDataSetFromHTSeq()

# ------------------------------------------
# Read in Prepared Annotation and Counts
countFiles = list.files(inDir, pattern="fb.txt$", full.names=TRUE)
basename(countFiles)
## [1] "treated1fb.txt" "treated2fb.txt" "treated3fb.txt" "untreated1fb.txt"
## [5] "untreated2fb.txt" "untreated3fb.txt" "untreated4fb.txt"
flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)
basename(flattenedFile)

# ------------------------------------------
# Standard Analysis Workflow

sample_info <- read.csv('data/sample_info.csv')
sample_info$CellType <- factor(sample_info$CellType)
sample_info$Condition <- factor(sample_info$Condition)

#Design taken for my DEseq script
design = ~ Condition + CellType

count_files <- c('data/HTSeq_Counts/_____',
                 'and etcetera')



# ------------------------------------------------------------------------------------
#  BELOW is just copied from Simone's DEXSeq Analysis (obviously needs modification before use here)
# ------------------------------------------------------------------------------------

# Create DEXSeq object
# Needed to remove quotation marks in .dexseq_count files (thanks to Arthur (@93355568) on bioconductor forums for solution)
dxd <- DEXSeqDataSetFromHTSeq(
  countfiles = count_files,  
  sampleData = sample_info, 
  design = ~ sample + exon + condition:exon,
  flattenedfile = 'data/gencode.v31.basic.annotation.gff'
)


# Normalization
dxd <- estimateSizeFactors(dxd)

# Dispersion Estimation (of the negative binomial distribution)
dxd <- estimateDispersions(dxd)
#plotDispEsts(dxd)


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
