library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ggupset)
library(clusterProfiler)
library(ggplot2)
library(stringr)
library(rtracklayer)
library(tidyverse)
library(ggrepel)
library(readxl)
library(biomaRt)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

ENHANCER_BED <- read.csv('data/CUTandTag/F5.hg38.enhancers.bed', sep='\t', header=F,
                         col.names = c('seqnames','start','end','name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts'))
ENHANCER_RANGES <- GRanges(seqnames=ENHANCER_BED$seqnames, ranges = IRanges(start = ENHANCER_BED$start, end = ENHANCER_BED$end))


# H3k27me3 peaks
h3k27me3_clustered <- read.csv('Results/CUTandTag/H3k27me3/csaw/csaw_sig_clustered_cpmFilt_k27me3_v5.bed',
                               sep='\t',col.names = c('seqnames','start','end','name','score','strand'),header=F)
h3k27me3_stats <- read.csv('Results/CUTandTag/H3k27me3/csaw/csaw_sig_clustered_cpmFilt_k27me3_STATS_v5.csv')
h3k27me3_grange <- GRanges(seqnames=h3k27me3_clustered$seqnames, ranges = IRanges(start = h3k27me3_clustered$start, end = h3k27me3_clustered$end),
                           direction=h3k27me3_stats$direction)
# H3k27ac peaks
h3k27ac_clustered <- read.csv('Results/CUTandTag/H3k27ac/csaw/csaw_sig_clustered_cpmFilt0.5_k27ac_v5.bed', 
                              sep='\t',col.names = c('seqnames','start','end','name','score','strand'),header=F)
h3k27ac_stats <- read.csv('Results/CUTandTag/H3k27ac/csaw/csaw_sig_clustered_cpmFilt_k27ac_STATS_v5.csv')
h3k27ac_grange <- GRanges(seqnames=h3k27ac_clustered$seqnames, ranges = IRanges(start = h3k27ac_clustered$start, end = h3k27ac_clustered$end),
                          direction=h3k27ac_stats$direction)
# H3k9me3 peaks
h3k9me3_clustered <- read.csv('Results/CUTandTag/H3k9me3/csaw/csaw_sig_clustered_cpmFilt_k9me3_v5.bed',
                               sep='\t',col.names = c('seqnames','start','end','name','score','strand'),header=F)
h3k9me3_stats <- read.csv('Results/CUTandTag/H3k9me3/csaw/csaw_sig_clustered_cpmFilt_k9me3_STATS_v5.csv')
h3k9me3_grange <- GRanges(seqnames=h3k9me3_clustered$seqnames, ranges = IRanges(start = h3k9me3_clustered$start, end = h3k9me3_clustered$end),
                           direction=h3k9me3_stats$direction)
h3k9me3_grange_merged100k <- read.csv('Results/CUTandTag/H3k9me3/csaw/csaw_k9me3_v6_merge50k_filt100k.bed',
                                      sep='\t',col.names = c('seqnames','start','end'),header=F)
h3k9me3_grange_merged100k_grange <- GRanges(seqnames=h3k9me3_grange_merged100k$seqnames, ranges = IRanges(start = h3k9me3_grange_merged100k$start, end = h3k9me3_grange_merged100k$end))
# ChIPseq peaks
prr12_chipseq <- read.csv('Results/CUTandTag/External/GSM8590159_THC0467WNHN_PRR12_ChIP2_hg38_Hamed_peaks.narrowPeak', sep='\t', header=F,
                         col.names = c('chrom','chromStart','chromEnd','name','score','strand','signalValue','pValue','qValue','peak'))
prr12_chipseq_filtered <- filter(prr12_chipseq, qValue>8)
prr12_grange <- GRanges(seqnames=prr12_chipseq$chrom, ranges = IRanges(start = prr12_chipseq$chromStart, end = prr12_chipseq$chromEnd))
prr12_filtered_grange <- GRanges(seqnames=prr12_chipseq_filtered$chrom, ranges = IRanges(start = prr12_chipseq_filtered$chromStart, end = prr12_chipseq_filtered$chromEnd))

res_NPC_KO_vs_WT <- read.csv('Results/DEGs/DEGs_NPC_KO_vs_WT.csv')
res_NPC_HET_vs_WT <- read.csv('Results/DEGs/DEGs_NPC_HET_vs_WT.csv')

NPC_KO_WT_sig_genes_up <- na.exclude(res_NPC_KO_vs_WT$X[res_NPC_KO_vs_WT$padj < 0.05 & res_NPC_KO_vs_WT$log2FoldChange>1])  # UP in KO
NPC_KO_WT_sig_genes_down <- na.exclude(res_NPC_KO_vs_WT$X[res_NPC_KO_vs_WT$padj < 0.05 & res_NPC_KO_vs_WT$log2FoldChange< -1]) # DOWN in KO


# Venn Diagram ------------------------------------------------------------
intersect(h3k27me3_grange,h3k27ac_grange,ignore.strand=T)
intersect(h3k27me3_grange,prr12_grange,ignore.strand=T)
intersect(h3k27ac_grange,prr12_grange,ignore.strand=T)


# ChIPseq Annotation ------------------------------------------------------
prr12_peakAnno <- annotatePeak(
  prr12_grange,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db",
  tssRegion = c(-3000, 3000)
)

plotAnnoPie(prr12_peakAnno)
plotDistToTSS(prr12_peakAnno, title="Distribution of Peaks relative to TSS")

upsetplot(prr12_peakAnno)

prr12_grange_promoter <- prr12_peakAnno@anno[prr12_peakAnno@detailGenomicAnnotation$Promoter,]
prr12_grange_enhancer <- GenomicRanges::intersect(prr12_grange, ENHANCER_RANGES, ignore.strand=TRUE)

prr12_enhancer_peakAnno <- annotatePeak(prr12_grange_enhancer,TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000))
prr12_enhancer_genes <- unique(prr12_enhancer_peakAnno@anno$SYMBOL)



# ChIPseq vs RNAseq -------------------------------------------------------



intersect(prr12_peakAnno@anno$SYMBOL, NPC_KO_WT_sig_genes_up)
intersect(prr12_peakAnno@anno$SYMBOL, NPC_KO_WT_sig_genes_down)

vennplot(list(prr12_chip=unique(prr12_peakAnno@anno$SYMBOL), RNAseq_UP=NPC_KO_WT_sig_genes_up, RNAseq_DOWN=NPC_KO_WT_sig_genes_down))



vennplot(list(h3k27me3=h3k27me3_grange, h3k27ac=h3k27ac_grange, h3k9me3=h3k9me3_grange, prr12=prr12_grange))

h3k27me3_grange_chip <- subsetByOverlaps(h3k27me3_grange, prr12_grange)
h3k27ac_grange_chip <- subsetByOverlaps(h3k27ac_grange, prr12_grange)
h3k9me3_grange_chip <- subsetByOverlaps(h3k9me3_grange, prr12_grange)
h3k9me3_grange_merged100k_grange_chip <- subsetByOverlaps(h3k9me3_grange_merged100k_grange, prr12_grange)

h3k27me3_grange_chip_peakAnno <- annotatePeak(h3k27me3_grange_chip,
                                              TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000))
h3k27ac_grange_chip_peakAnno <- annotatePeak(h3k27ac_grange_chip,
                                              TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000))
h3k9me3_grange_chip_peakAnno <- annotatePeak(h3k9me3_grange_chip,
                                              TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000))
h3k9me3_grange_merged100k_grange_chip_peakAnno <- annotatePeak(h3k9me3_grange_merged100k_grange_chip,
                                              TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000))

h3k27me3_grange_chip_genes <- unique(h3k27me3_grange_chip_peakAnno@anno$SYMBOL)
h3k27ac_grange_chip_genes <- unique(h3k27ac_grange_chip_peakAnno@anno$SYMBOL)
h3k9me3_grange_chip_genes <- unique(h3k9me3_grange_chip_peakAnno@anno$SYMBOL)
h3k9me3_grange_merged100k_grange_chip_genes <- unique(h3k9me3_grange_merged100k_grange_chip_peakAnno@anno$SYMBOL)

intersect(NPC_KO_WT_sig_genes_up, h3k27me3_grange_chip_genes)
intersect(NPC_KO_WT_sig_genes_down, h3k27me3_grange_chip_genes)
intersect(NPC_KO_WT_sig_genes_up, h3k27ac_grange_chip_genes)
intersect(NPC_KO_WT_sig_genes_down, h3k27ac_grange_chip_genes)
intersect(NPC_KO_WT_sig_genes_up, h3k9me3_grange_chip_genes)
intersect(NPC_KO_WT_sig_genes_down, h3k9me3_grange_chip_genes)
intersect(NPC_KO_WT_sig_genes_up, h3k9me3_grange_merged100k_grange_chip_genes)
intersect(NPC_KO_WT_sig_genes_down, h3k9me3_grange_merged100k_grange_chip_genes)



# Compare Manual H3k9me3 region genes -------------------------------------
higherKO <- c(
  "AGBL4",
  "TCHH", "TCHHL1", "RPTN", "HRNR", "CCDST",
  "RP11-157E21.1", "LOC101927798",
  "FAM135B", "COL22A1", "KCNK9", "TRAPPC9",
  "OR1J2", "OR1N1", "OR1N2", "OR1L8", "OR1Q1", "OR1B1", "OR1L3", "OR1L4", "OR1L6",
  "CACNA1B",
  "CYP4F22", "CYP4F8", "CYP4F3", "CYP4F12", "OR10H2", "CYP4F24P", "OR10H5",
  "UCA1", "CLEC4OP", "CYP4F2", "CYP4F11", "OR10H4", "TPM4", "PGLYRP2", "RASAL3", "WIZ",
  "LIPE-AS1", "PSG3", "PSG8",
  "SIGLEC8", "CEACAM18", "SIGLEC12", "SIGLEC6", "ZNF175", "SIGLEC5", "SIGLEC14",
  "DEFB115", "DEFB116", "DEFB118", "DEFB119", "DEFB121", "DEFB122", "DEFB123", "DEFB124",
  "REM1", "HM13",
  "MC3R"
)

lowerKO <- c(
  "AJAP1",
  "SLURP1", "LYPD2", "LYNX1", "LY6D", "GML",
  "ZNF132", "ZNF324B", "ZNF324", "ZNF446", "SLC27A5", "ZBTB45", "TRIM28",
  "BCOR", "RP11-320G24.1"
)


# PRR12 ChIPseq Filtering -------------------------------------------------
library(readxl)

prr12_chip_excel <- read_xlsx('Results/CUTandTag/External/GSM8590159_THC0467WNHN_PRR12_ChIP2_hg38_Hamed_peaks.narrowPeak.xlsx')


prr12_chip_excel_filt <- filter(prr12_chip_excel, qValue>8)

write.table(prr12_chip_excel_filt[,1:3], 'Results/CUTandTag/External/PRR12_ChIP2_filtered.narrowPeak', quote = F, sep='\t', col.names = F, row.names = F)


# Histone Mods W/ Direction x PRR12 ChIP -----------------------------------------------

h3k27me3_grange_up <- h3k27me3_grange[h3k27me3_stats$direction=='up']
h3k27me3_grange_down <- h3k27me3_grange[h3k27me3_stats$direction=='down']
h3k27ac_grange_up <- h3k27ac_grange[h3k27ac_stats$direction=='up']
h3k27ac_grange_down <- h3k27ac_grange[h3k27ac_stats$direction=='down']
h3k9me3_grange_up <- h3k9me3_grange[h3k9me3_stats$direction=='up']
h3k9me3_grange_down <- h3k9me3_grange[h3k9me3_stats$direction=='down']

h3k27me3_grange_up_chip <- GenomicRanges::intersect(h3k27me3_grange_up, prr12_grange)
h3k27me3_grange_down_chip <- GenomicRanges::intersect(h3k27me3_grange_down, prr12_grange)
h3k27ac_grange_up_chip <- GenomicRanges::intersect(h3k27ac_grange_up, prr12_grange)
h3k27ac_grange_down_chip <- GenomicRanges::intersect(h3k27ac_grange_down, prr12_grange)
h3k9me3_grange_up_chip <- GenomicRanges::intersect(h3k9me3_grange_up, prr12_grange)
h3k9me3_grange_down_chip <- GenomicRanges::intersect(h3k9me3_grange_down, prr12_grange)



h3k27me3_grange_up_peakAnno <- annotatePeak(h3k27me3_grange_up_chip,TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000))
h3k27me3_grange_down_peakAnno <- annotatePeak(h3k27me3_grange_down_chip,TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000))
h3k27ac_grange_up_peakAnno <- annotatePeak(h3k27ac_grange_up_chip,TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000))
h3k27ac_grange_down_peakAnno <- annotatePeak(h3k27ac_grange_down_chip,TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000))
h3k9me3_grange_up_peakAnno <- annotatePeak(h3k9me3_grange_up_chip,TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000))
h3k9me3_grange_down_peakAnno <- annotatePeak(h3k9me3_grange_down_chip,TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000))



# PRR12 x H3K27ac ------------------------------------------------------------

# H3k27ac ChIP peaks
h3k27ac_chip <- read.csv('Results/CUTandTag/External/ENCODE_or_ChIPAtlas/h3k27ac_ENCFF225FPN.bed', sep = '\t', header=F,
                          col.names = c('chrom','chromStart','chromEnd','name','score','strand','signalValue','pValue','qValue','peak'))
h3k27ac_chip_grange <- GRanges(seqnames=h3k27ac_chip$chrom, ranges = IRanges(start = h3k27ac_chip$chromStart, end = h3k27ac_chip$chromEnd))



prr12_h3k27acChip_ranges <- GenomicRanges::intersect(prr12_filtered_grange, h3k27ac_chip_grange)

prr12_h3k27ac_peakAnno <- annotatePeak(
  prr12_h3k27acChip_ranges,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db",
  tssRegion = c(-3000, 3000)
)

plotAnnoPie(prr12_h3k27ac_peakAnno, legend.position = 'ta')


prr12_h3k27acChip_ranges_enhancer <- GenomicRanges::intersect(prr12_h3k27acChip_ranges, ENHANCER_RANGES, ignore.strand=T)


# custom piechart with enhancers
anno_df <- as.data.frame(prr12_h3k27ac_peakAnno)
enh_hits <- overlapsAny(
  prr12_h3k27acChip_ranges,
  ENHANCER_RANGES,
  ignore.strand = TRUE
)
anno_df$custom_annotation <- anno_df$annotation
anno_df$custom_annotation[enh_hits] <- "Enhancer"
anno_df$custom_annotation <- gsub("\\s*\\(.*\\)", "", anno_df$custom_annotation)

pie_df <- anno_df %>%
  count(custom_annotation) %>%
  arrange(desc(n)) %>%
  mutate(
    fraction = n / sum(n),
    ymax = cumsum(fraction),
    ymin = lag(ymax, default = 0),
    ymid = (ymax + ymin) / 2,
    label = paste0(
      custom_annotation, ' (',
      round(fraction * 100, 1),
      "%)"
    )
  )
ggplot(pie_df, aes(x = 1, y = fraction, fill = custom_annotation)) +
  geom_col(color = "white", width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.2) +   # extra room for labels
  theme_void() +
  theme(
    legend.position = "none"
  ) +
  geom_label_repel(
    aes(
      x = 1.3,
      y = ymid,
      label = label
    ),
    size = 4,
    nudge_x = 0.6,
    show.legend = FALSE,
    segment.color = "grey50",
    segment.size = 0.7,
    direction = "y",
    hjust = 0
  ) 

promoter_hits <- grepl("Promoter", anno_df$annotation)
prr12_h3k27acChip_promOrEnhancer <- prr12_h3k27acChip_ranges[
  overlapsAny(prr12_h3k27acChip_ranges, ENHANCER_RANGES,ignore.strand=T) | promoter_hits]
prr12_h3k27acChip_Enhancer <- prr12_h3k27acChip_ranges[
  overlapsAny(prr12_h3k27acChip_ranges, ENHANCER_RANGES,ignore.strand=T)]

prr12_h3k27acChip_promOrEnhancer_peakAnno <- annotatePeak(
  prr12_h3k27acChip_promOrEnhancer,TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000))
prr12_h3k27acChip_Enhancer_peakAnno <- annotatePeak(
  prr12_h3k27acChip_Enhancer,TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000))

prr12_h3k27acChip_promOrEnhancer_genes <- unique(as.data.frame(prr12_h3k27acChip_promOrEnhancer_peakAnno)$geneId)
prr12_h3k27acChip_Enhancer_genes <- unique(as.data.frame(prr12_h3k27acChip_Enhancer_peakAnno)$geneId)

# Prep Background
background_genelist <- unique(
  as.data.frame(prr12_h3k27ac_peakAnno)$geneId
)

# GO enrichment
ego <- enrichGO(
  gene          = prr12_h3k27acChip_Enhancer_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",      # Biological Process
  pAdjustMethod = "BH",
  universe = background_genelist,
  qvalueCutoff  = 0.05,
  readable = T
)

barplot(ego, showCategory=18)

dotplot(ego, showCategory = 15) + ggtitle(NULL) + 
  theme(text=element_text(size=24), axis.text.x = element_text(size=24), axis.text.y = element_text(size=24), axis.title.x = element_text(size=24)) + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))


## Extract Term Genes
ego@result
  filter(ego@result, grepl('eye',ego@result$Description))


GO_TERM_OF_INTEREST <- 'eye development'
go_genes <- strsplit(ego@result$geneID[ego@result$Description==GO_TERM_OF_INTEREST], split = '/')[[1]]
writeClipboard(go_genes)



# -------------------------------------------------------------------------
# Cross GO Terms ----------------------------------------------------------

goterm_genes_xl <- read_xlsx('Results/GOTerm_Genes_CUTandTag.xlsx', sheet = 'ChIP PRR12 x H3K27AC x PromEnha')

goterm_genes_unique <- as.character(na.exclude(unique(as.vector(as.matrix(goterm_genes_xl)))))

prr12chip_h3k27acDiff <- prr12_filtered_grange[overlapsAny(prr12_filtered_grange, h3k27ac_grange)]

prr12chip_h3k27acDiff_peakAnno <- annotatePeak(prr12chip_h3k27acDiff,TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000), addFlankGeneInfo = T)

# Extract and convert flank gene IDs
flank_genes <- c(na.omit(unique(unlist(lapply(prr12chip_h3k27acDiff_peakAnno@anno$flank_geneIds, strsplit, split=';')))))
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)
flank_genes_symbols <- getBM(
  attributes = c(
    "entrezgene_id",
    "hgnc_symbol",
    "ensembl_gene_id",
    "description"
  ),
  filters = "entrezgene_id",
  values = flank_genes,
  mart = mart
)

# Continue with both 'main' genes and the flank genes
prr12chip_h3k27ac_genes <- unique(c(prr12chip_h3k27acDiff_peakAnno@anno$SYMBOL, flank_genes_symbols$hgnc_symbol))

prr12chip_h3k27ac_GO_genes <- intersect(prr12chip_h3k27ac_genes, goterm_genes_unique)

prr12chip_h3k27ac_GO_genes[prr12chip_h3k27ac_GO_genes %in% c(NPC_KO_WT_sig_genes_up, NPC_KO_WT_sig_genes_down)]














