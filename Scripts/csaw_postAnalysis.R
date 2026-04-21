library(csaw)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ggupset)
library(clusterProfiler)
library(ggplot2)
library(stringr)
library(rtracklayer)

ENHANCER_BED <- read.csv('data/CUTandTag/F5.hg38.enhancers.bed', sep='\t', header=F,
                         col.names = c('seqnames','start','end','name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts'))
ENHANCER_RANGES <- GRanges(seqnames=ENHANCER_BED$seqnames, ranges = IRanges(start = ENHANCER_BED$start, end = ENHANCER_BED$end))

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene



HISTONE_MODIFICATION <- 'H3k27me3'
HISTONE_MODIFICATION <- 'H3k27ac'
HISTONE_MODIFICATION <- 'H3k9me3'

if(HISTONE_MODIFICATION == 'H3k27me3'){
  bed_merged <- read.csv('Results/CUTandTag/H3k27me3/csaw/csaw_sig_merged_k27me3_v5.bed', sep='\t',col.names = c('seqnames','start','end','name','score','strand'),header=F)
  bed_clustered <- read.csv('Results/CUTandTag/H3k27me3/csaw/csaw_sig_clustered_cpmFilt_k27me3_v5.bed',
                            sep='\t',col.names = c('seqnames','start','end','name','score','strand'),header=F)
  clustered_stats <- read.csv('Results/CUTandTag/H3k27me3/csaw/csaw_sig_clustered_cpmFilt_k27me3_STATS_v5.csv')
} else if(HISTONE_MODIFICATION == 'H3k27ac'){
  bed_merged <- read.csv('Results/CUTandTag/H3k27ac/csaw/csaw_sig_merged_k27ac_v5.bed', sep='\t',col.names = c('seqnames','start','end','name','score','strand'),header=F)
  bed_clustered <- read.csv('Results/CUTandTag/H3k27ac/csaw/csaw_sig_clustered_cpmFilt_k27ac_v5.bed', 
                            sep='\t',col.names = c('seqnames','start','end','name','score','strand'),header=F)
  clustered_stats <- read.csv('Results/CUTandTag/H3k27ac/csaw/csaw_sig_clustered_cpmFilt_k27ac_STATS_v5.csv')
} else if(HISTONE_MODIFICATION == 'H3k9me3'){
  bed_merged <- read.csv('Results/CUTandTag/H3k9me3/csaw/csaw_sig_merged_k9me3_v5.bed', sep='\t',col.names = c('seqnames','start','end','name','score','strand'),header=F)
  bed_clustered <- read.csv('Results/CUTandTag/H3k9me3/csaw/csaw_sig_clustered_cpmFilt_k9me3_v5.bed',
                            sep='\t',col.names = c('seqnames','start','end','name','score','strand'),header=F)
  clustered_stats <- read.csv('Results/CUTandTag/H3k9me3/csaw/csaw_sig_clustered_cpmFilt_k9me3_STATS_v5.csv')
}

grange_merged <- GRanges(seqnames=bed_merged$seqnames, ranges = IRanges(start = bed_merged$start, end = bed_merged$end))
grange_clustered <- GRanges(seqnames=bed_clustered$seqnames, ranges = IRanges(start = bed_clustered$start, end = bed_clustered$end))

grange_selected <- grange_clustered





peakAnno <- annotatePeak(
  grange_selected,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db",
  tssRegion = c(-3000, 3000)
)

plotAnnoPie(peakAnno)
plotDistToTSS(peakAnno, title="Distribution of Peaks relative to TSS")

upsetplot(peakAnno)

# Load and annotate background
if(HISTONE_MODIFICATION == 'H3k27me3'){
  
} else if(HISTONE_MODIFICATION == 'H3k27ac'){
  background_bed <- read.csv('Results/CUTandTag/H3k27ac/epic2_multiinter_supp5.bed', sep='\t',col.names = c('seqnames','start','end'))
} else if(HISTONE_MODIFICATION == 'H3k9me3'){
  
}

background_granges <- GRanges(seqnames=background_bed$seqnames, ranges = IRanges(start = background_bed$start,end = background_bed$end))

background_peakAnno <- annotatePeak(
  background_granges,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db",
  tssRegion = c(-3000, 3000)
)

background_genelist <- unique(as.data.frame(background_peakAnno)$geneId)

# --------  GO Analysis taken from online
# Get gene IDs from annotation
gene_list <- unique(as.data.frame(peakAnno)$geneId)

# GO enrichment
ego <- enrichGO(
  gene          = gene_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",      # Biological Process
  pAdjustMethod = "BH",
#  universe = background_genelist,
  qvalueCutoff  = 0.05,
  readable = T
)

# KEGG pathway enrichment
#kk <- enrichKEGG(gene = gene_list, organism = "hsa")

# Visualize top GO terms
barplot(ego, showCategory=20, title="GO Enrichment")
#dotplot(kk, showCategory=20, title="KEGG Enrichment")

dotplot(ego, showCategory = 20) + ggtitle(NULL) + 
  theme(text=element_text(size=24), axis.text.x = element_text(size=24), axis.text.y = element_text(size=24), axis.title.x = element_text(size=24)) + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))

# ----- For multiple conditions
# List of annotated peaks for two conditions
peakAnno_list <- list(ConditionA=peakAnnoA, ConditionB=peakAnnoB)

# Compare feature distribution
plotAnnoBar(peakAnno_list)




# -------------------------------------------------------------------------------
####### Compare Across Histone Modifications ------------------------------------

h3k27me3_clustered <- read.csv('Results/CUTandTag/H3k27me3/csaw/csaw_sig_clustered_cpmFilt_k27me3_v5.bed',
                          sep='\t',col.names = c('seqnames','start','end','name','score','strand'),header=F)
h3k27me3_stats <- read.csv('Results/CUTandTag/H3k27me3/csaw/csaw_sig_clustered_cpmFilt_k27me3_STATS_v5.csv')
h3k27me3_grange <- GRanges(seqnames=h3k27me3_clustered$seqnames, ranges = IRanges(start = h3k27me3_clustered$start, end = h3k27me3_clustered$end),
                           direction=h3k27me3_stats$direction)
overlap_hits <- findOverlaps(h3k27me3_grange, ENHANCER_RANGES)
h3k27me3_grange$enhancer <- FALSE
h3k27me3_grange$enhancer[queryHits(overlap_hits)] <- TRUE
h3k27me3_peakAnno <- annotatePeak(h3k27me3_grange,TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000),level = 'transcript')
genelist_h3k27me3 <- unique(h3k27me3_peakAnno@anno$SYMBOL)
genelist_h3k27me3_up <- unique(h3k27me3_peakAnno@anno$SYMBOL[h3k27me3_peakAnno@anno$direction=='up'])
genelist_h3k27me3_down <- unique(h3k27me3_peakAnno@anno$SYMBOL[h3k27me3_peakAnno@anno$direction=='down'])

h3k27ac_clustered <- read.csv('Results/CUTandTag/H3k27ac/csaw/csaw_sig_clustered_cpmFilt_k27ac_v5.bed', 
                          sep='\t',col.names = c('seqnames','start','end','name','score','strand'),header=F)
h3k27ac_stats <- read.csv('Results/CUTandTag/H3k27ac/csaw/csaw_sig_clustered_cpmFilt_k27ac_STATS_v5.csv')
h3k27ac_grange <- GRanges(seqnames=h3k27ac_clustered$seqnames, ranges = IRanges(start = h3k27ac_clustered$start, end = h3k27ac_clustered$end),
                          direction=h3k27ac_stats$direction)
overlap_hits <- findOverlaps(h3k27ac_grange, ENHANCER_RANGES)
h3k27ac_grange$enhancer <- FALSE
h3k27ac_grange$enhancer[queryHits(overlap_hits)] <- TRUE
h3k27ac_peakAnno <- annotatePeak(h3k27ac_grange,TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000),level = 'transcript')
genelist_h3k27ac <- unique(h3k27ac_peakAnno@anno$SYMBOL)
genelist_h3k27ac_up <- unique(h3k27ac_peakAnno@anno$SYMBOL[h3k27ac_peakAnno@anno$direction=='up'])
genelist_h3k27ac_down <- unique(h3k27ac_peakAnno@anno$SYMBOL[h3k27ac_peakAnno@anno$direction=='down'])

h3k9me3_clustered <- read.csv('Results/CUTandTag/H3k9me3/csaw/csaw_sig_clustered_cpmFilt_k9me3_v5.bed',
                          sep='\t',col.names = c('seqnames','start','end','name','score','strand'),header=F)
h3k9me3_stats <- read.csv('Results/CUTandTag/H3k9me3/csaw/csaw_sig_clustered_cpmFilt_k9me3_STATS_v5.csv')
h3k9me3_grange <- GRanges(seqnames=h3k9me3_clustered$seqnames, ranges = IRanges(start = h3k9me3_clustered$start, end = h3k9me3_clustered$end),
                          direction=h3k9me3_stats$direction)
overlap_hits <- findOverlaps(h3k9me3_grange, ENHANCER_RANGES)
h3k9me3_grange$enhancer <- FALSE
h3k9me3_grange$enhancer[queryHits(overlap_hits)] <- TRUE
h3k9me3_peakAnno <- annotatePeak(h3k9me3_grange,TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000),level = 'transcript')
genelist_h3k9me3 <- unique(h3k9me3_peakAnno@anno$SYMBOL)
genelist_h3k9me3_up <- unique(h3k9me3_peakAnno@anno$SYMBOL[h3k9me3_peakAnno@anno$direction=='up'])
genelist_h3k9me3_down <- unique(h3k9me3_peakAnno@anno$SYMBOL[h3k9me3_peakAnno@anno$direction=='down'])
h3k9me3_background_ranges <- readRDS('Results/CUTandTag/H3k9me3/csaw/universe_ranges_WTvKO_v6.RData')
h3k9me3_BackgroundAnno <- annotatePeak(h3k9me3_background_ranges,TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000),level = 'transcript')
h3k9me3_background_genes <- unique(h3k9me3_BackgroundAnno@anno$SYMBOL)


vennplot(list(k27me3_up=genelist_h3k27me3_up, k27me3_down=genelist_h3k27me3_down,
              ky27ac_up=genelist_h3k27ac_up, k27ac_down=genelist_h3k27ac_down))


### H3k27ac peaks that have a gene shared with H3k27me3
# H3k27me3 upreg and promoter, H3k27ac downreg and promoter/enhancer
h3k27ac_anno_cross_subset <- h3k27ac_peakAnno@anno[h3k27ac_peakAnno@anno$SYMBOL %in% genelist_h3k27me3_up,]
h3k27ac_peakAnno_subet <- h3k27ac_peakAnno
h3k27ac_peakAnno_subet@anno <- h3k27ac_anno_cross_subset
h3k27ac_subset_grange <- GRanges(seqnames=h3k27ac_anno_cross_subset@seqnames, ranges = h3k27ac_anno_cross_subset@ranges,
                           direction=h3k27ac_anno_cross_subset$direction)
overlap_hits <- findOverlaps(h3k27ac_subset_grange, ENHANCER_RANGES)
h3k27ac_subset_grange$enhancer <- FALSE
h3k27ac_subset_grange$enhancer[queryHits(overlap_hits)] <- TRUE
h3k27ac_subset_peakAnno <- annotatePeak(h3k27ac_subset_grange,TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000))

h3k27me3_upProm_h3k27ac_downEnhancer <- h3k27ac_subset_peakAnno@anno$SYMBOL[h3k27ac_subset_peakAnno@anno$enhancer==TRUE & h3k27ac_subset_peakAnno@anno$direction=='down']
h3k27me3_upProm_h3k27ac_downProm <- intersect(h3k27ac_subset_peakAnno@anno$SYMBOL[grep('Promoter',h3k27ac_subset_peakAnno@anno$annotation)],
          h3k27ac_subset_peakAnno@anno$SYMBOL[h3k27ac_subset_peakAnno@anno$direction=='down'])
# H3k27me3 downreg and promoter, H3k27ac upreg and promoter/enhancer
h3k27ac_anno_cross_subset <- h3k27ac_peakAnno@anno[h3k27ac_peakAnno@anno$SYMBOL %in% genelist_h3k27me3_down,]
h3k27ac_peakAnno_subet <- h3k27ac_peakAnno
h3k27ac_peakAnno_subet@anno <- h3k27ac_anno_cross_subset
h3k27ac_subset_grange <- GRanges(seqnames=h3k27ac_anno_cross_subset@seqnames, ranges = h3k27ac_anno_cross_subset@ranges,
                                 direction=h3k27ac_anno_cross_subset$direction)
overlap_hits <- findOverlaps(h3k27ac_subset_grange, ENHANCER_RANGES)
h3k27ac_subset_grange$enhancer <- FALSE
h3k27ac_subset_grange$enhancer[queryHits(overlap_hits)] <- TRUE
h3k27ac_subset_peakAnno <- annotatePeak(h3k27ac_subset_grange,TxDb = txdb,annoDb = "org.Hs.eg.db",tssRegion = c(-3000, 3000))

h3k27me3_downProm_h3k27ac_upEnhancer <- h3k27ac_subset_peakAnno@anno$SYMBOL[h3k27ac_subset_peakAnno@anno$enhancer==TRUE & h3k27ac_subset_peakAnno@anno$direction=='up']
h3k27me3_downProm_h3k27ac_upProm <- intersect(h3k27ac_subset_peakAnno@anno$SYMBOL[grep('Promoter',h3k27ac_subset_peakAnno@anno$annotation)],
                                                  h3k27ac_subset_peakAnno@anno$SYMBOL[h3k27ac_subset_peakAnno@anno$direction=='up'])


## Cross with RNAseq
res_NPC_KO_vs_WT <- read.csv('Results/DEGs/DEGs_NPC_KO_vs_WT.csv')
res_NPC_HET_vs_WT <- read.csv('Results/DEGs/DEGs_NPC_HET_vs_WT.csv')

NPC_KO_WT_sig_genes_up <- na.exclude(res_NPC_KO_vs_WT$X[res_NPC_KO_vs_WT$padj < 0.05 & res_NPC_KO_vs_WT$log2FoldChange>0])  # UP in KO
NPC_KO_WT_sig_genes_down <- na.exclude(res_NPC_KO_vs_WT$X[res_NPC_KO_vs_WT$padj < 0.05 & res_NPC_KO_vs_WT$log2FoldChange<0]) # DOWN in KO


intersect(h3k27me3_upProm_h3k27ac_downEnhancer, NPC_KO_WT_sig_genes_down)
intersect(h3k27me3_upProm_h3k27ac_downProm, NPC_KO_WT_sig_genes_down)

intersect(h3k27me3_downProm_h3k27ac_upEnhancer, NPC_KO_WT_sig_genes_up)
intersect(h3k27me3_downProm_h3k27ac_upProm, NPC_KO_WT_sig_genes_up)

### H3k27me3 peaks crossed with RNAseq


### H3k9me3 peaks crossed with RNAseq
vennplot(list(h3k9me3_peak=h3k9me3_peakAnno@anno$SYMBOL, upreg=NPC_KO_WT_sig_genes_up, downreg=NPC_KO_WT_sig_genes_down))
vennplot(list(h3k9me3_peak_upKO=h3k9me3_peakAnno@anno$SYMBOL[h3k9me3_peakAnno@anno$direction=='up'],
              h3k9me3_peak_downKO=h3k9me3_peakAnno@anno$SYMBOL[h3k9me3_peakAnno@anno$direction=='down'],
              upreg=NPC_KO_WT_sig_genes_up, downreg=NPC_KO_WT_sig_genes_down))

h3k9me3_peak_and_NPC_sig_up <- intersect(h3k9me3_peakAnno@anno$SYMBOL, NPC_KO_WT_sig_genes_up)
h3k9me3_peak_and_NPC_sig_down <- intersect(h3k9me3_peakAnno@anno$SYMBOL, NPC_KO_WT_sig_genes_down)


h3k9me3_peak_upKO_and_NPC_sig_down <- intersect(h3k9me3_peakAnno@anno$SYMBOL[h3k9me3_peakAnno@anno$direction=='up'], NPC_KO_WT_sig_genes_down)
h3k9me3_peak_downKO_and_NPC_sig_up <- intersect(h3k9me3_peakAnno@anno$SYMBOL[h3k9me3_peakAnno@anno$direction=='down'], NPC_KO_WT_sig_genes_up)


ego <- enrichGO(
  gene          = h3k9me3_peak_upKO_and_NPC_sig_down,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",      # Biological Process
  pAdjustMethod = "BH",
  universe = intersect(res_NPC_KO_vs_WT$X[res_NPC_KO_vs_WT$baseMean>5], h3k9me3_background_genes),
  qvalueCutoff  = 0.05,
  readable = T
)
barplot(ego, showCategory=20, title="GO Enrichment")


#### Extract H3k27me3 promoter peaks / H3k27ac promoter+enhancer peaks

h3k27me3_promoters <- h3k27me3_grange[h3k27me3_peakAnno@detailGenomicAnnotation$Promoter,]

export(h3k27me3_promoters, 'Results/CUTandTag/H3k27me3/csaw/csaw_sig_clustered_cpmFilt2_k27me3_PROMOTERS_v5.bed')




h3k27ac_promoters_enhancers <- h3k27ac_grange[h3k27ac_peakAnno@detailGenomicAnnotation$Promoter | h3k27ac_grange$enhancer,]
export(h3k27ac_promoters_enhancers, 'Results/CUTandTag/H3k27ac/csaw/csaw_sig_clustered_cpmFilt0.5_k27ac_PROMandENHA_v5.bed')


### Cross H3k9me3 and RNAseq with homeodomain genes
homeodomain_factors_db <- read.csv('/Users/tbehr/Desktop/H14CORE_motifs.tsv', sep='\t')
homeodomain_genes <- unique(homeodomain_factors_db$Gene..human.)

vennplot(list(h3k9me3_peak_upKO=h3k9me3_peakAnno@anno$SYMBOL[h3k9me3_peakAnno@anno$direction=='up'],
              h3k9me3_peak_downKO=h3k9me3_peakAnno@anno$SYMBOL[h3k9me3_peakAnno@anno$direction=='down'],
              homeodomain_genes=homeodomain_genes))

vennplot(list(NPC_KO_WT_sig_genes_up=NPC_KO_WT_sig_genes_up,
              NPC_KO_WT_sig_genes_down=NPC_KO_WT_sig_genes_down,
              homeodomain_genes=homeodomain_genes))





genes_custom <- c("AJAP1",
           "AGBL4",
           "TCHH", "TCHHL1", "RPTN", "HRNR", "CCDST",
           "RP11-157E21.1", "LOC101927798",
           "FAM135B", "COL22A1", "KCNK9", "TRAPPC9",
           "SLURP1", "LYPD2", "LYNX1", "LY6D", "GML",
           "OR1J2", "OR1N1", "OR1N2", "OR1L8", "OR1Q1", "OR1B1", "OR1L3", "OR1L4", "OR1L6",
           "CACNA1B",
           "CYP4F22", "CYP4F8", "CYP4F3", "CYP4F12", "OR10H2", "CYP4F24P", "OR10H5", "UCA1", "CLEC4OP",
           "CYP4F2", "CYP4F11", "OR10H4", "TPM4", "PGLYRP2", "RASAL3", "WIZ",
           "LIPE-AS1", "PSG3", "PSG8",
           "SIGLEC8", "CEACAM18", "SIGLEC12", "SIGLEC6", "ZNF175", "SIGLEC5", "SIGLEC14",
           "ZNF132", "ZNF324B", "ZNF324", "ZNF446", "SLC27A5", "ZBTB45", "TRIM28",
           "DEFB115", "DEFB116", "DEFB118", "DEFB119", "DEFB121", "DEFB122", "DEFB123", "DEFB124", "REM1", "HM13",
           "MC3R",
           "BCOR", "RP11-320G24.1")

