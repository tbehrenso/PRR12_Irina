library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggbreak)

CODEP_TABLE <- read.csv("data/External/PRR12's Top 100 Codependencies for CRISPR (DepMap Public 26Q1+Score, Chronos).csv")
CODEP_TABLE$rank <- as.numeric(rownames(CODEP_TABLE))
GENES_TO_LABEL <- c('KAT2A','HDAC8','EP300','EHMT1','EHMT2','BRD9','SMARCD1','TADA1','USP22')
GENES_TO_LABEL_SUPPLEMENTARY <- c("TADA1","TADA2B","TADA3","SUPT7L","SUPT20H","USP22","SGF29","TAF5L","TAF6L")

label_data <- data.frame(Gene=GENES_TO_LABEL,
                         Correlation=CODEP_TABLE$Correlation[match(GENES_TO_LABEL, CODEP_TABLE$Gene)],
                         label = GENES_TO_LABEL,
                         label_color = rep('#ff0000',length(GENES_TO_LABEL)),
                         rank = CODEP_TABLE$rank[match(GENES_TO_LABEL, CODEP_TABLE$Gene)])

ggplot(CODEP_TABLE, aes(x=rank,y=Correlation)) + geom_point() + theme_minimal() + theme(legend.position='none',panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0.16, linetype='dashed') +
  geom_hline(yintercept = -0.16, linetype='dashed') +
  xlab('Rank') + ylab('Co-dependency Correlation') +
  geom_text_repel(
    data = label_data,
    aes(label = label, color = label_color),
    size = 6,
    show.legend = FALSE,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    fontface = "bold",
    nudge_y = c(-0.08,0.08,0.08,-0.08,0.08,0.08,0.08,0.08,-0.08), nudge_x = 5) + 
  ylim(c(-0.2,0.6)) +
  geom_hline(yintercept = 0, linetype='solid', linewidth=1.5) +
  geom_vline(xintercept = 0, linetype='solid', linewidth=1.5)
