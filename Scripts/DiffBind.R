library(DiffBind)
library(plyranges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(org.Hs.eg.db)

samplesheet <- read.csv('data/samplesheet_diffbind.csv')

#db <- dba(sampleSheet = samplesheet)
#db <- dba.count(db, summits = 250)

#db <- readRDS('data/db_diff.RData')
#db <- readRDS('data/db_normalized_diff.RData')
#db <- readRDS('data/db_reduced_diff.RData')
db <- readRDS('data/db_reduced_normalized_diff.RData')


dba.plotPCA(db,  attributes=DBA_CONDITION, label=DBA_ID)

# db <- dba.normalize(db)

# # define the contrast
# db <- dba.contrast(
#   db,
#   categories = DBA_CONDITION,
#   minMembers = 2
# )
# 
# # differential analysis
# db <- dba.analyze(
#   db,
#   method = DBA_DESEQ2
# )

dba.show(db, bContrasts=T)

dba.plotMA(db, method=DBA_DESEQ2, contrast = 2)
dba.plotMA(db, bXY=TRUE, contrast = 2)
dba.plotVolcano(db, contrast = 2)

pvals <- dba.plotBox(db, contrast = 2)


res_deseq <- dba.report(db, method = DBA_DESEQ2, contrast = 2, th = 1)

res_deseq_filt <- res_deseq %>%
  filter(Fold > 0, FDR < 0.05)

res_deseq_df <- as.data.frame(res_deseq)
hist(res_deseq_df$Fold)
res_deseq_df[res_deseq_df$FDR<0.05,]



# ------------------
# Annotate with Genes

#h38_gtf <- makeTxDbFromGFF('/Users/tbehr/Desktop/SanRaffaele/Projects/REFERENCE/Homo_sapiens.GRCh38.115.gtf')
h38_gtf <- import('/Users/tbehr/Desktop/SanRaffaele/Projects/REFERENCE/Homo_sapiens.GRCh38.115.gtf')
#h38_gtf <- import('/Users/tbehr/Desktop/SanRaffaele/Projects/REFERENCE/NCBI_ref_hg38_GCF_000001405.40.gtf')

seqlevels(h38_gtf) <- paste0("chr", seqlevels(h38_gtf))

annotation_matches <- findOverlaps(res_deseq, h38_gtf, maxgap = 5000)

res_deseq$gene_id <- NA_character_
res_deseq$gene_id[queryHits(annotation_matches)] <- h38_gtf$gene_name[subjectHits(annotation_matches)]

# Annotate with genomic feature
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peak_anno <- annotatePeak(
  res_deseq,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)
peak_anno_granges <- as.GRanges(peak_anno)

peak_anno_granges$anno_simple <- dplyr::case_when(
  grepl("^Promoter", peak_anno_granges$annotation) ~ "Promoter",
  grepl("^Exon", peak_anno_granges$annotation) ~ "Exon",
  grepl("^Intron", peak_anno_granges$annotation) ~ "Intron",
  grepl("5' UTR", peak_anno_granges$annotation) ~ "5' UTR",
  grepl("3' UTR", peak_anno_granges$annotation) ~ "3' UTR",
  grepl("Downstream", peak_anno_granges$annotation) ~ "Downstream",
  grepl("Distal Intergenic", peak_anno_granges$annotation) ~ "Distal Intergenic",
  TRUE ~ "Other"
)

# correct very distant TSS annotations
peak_anno_granges$SYMBOL[abs(peak_anno_granges$distanceToTSS) > 20000] <- NA

res_deseq <- peak_anno_granges

if(PROMOTER_PEAKS == TRUE){
  res_deseq <- res_deseq[res_deseq$anno_simple == 'Promoter',]
}


diff_binded_genes <- unique(na.exclude(res_deseq$SYMBOL[res_deseq$FDR < 0.05]))
#diff_binded_genes <- unique(na.exclude(res_deseq$SYMBOL[res_deseq$FDR < 0.05 & abs(res_deseq$Fold) > 1]))
background_genes <- unique(na.exclude(res_deseq$SYMBOL))

# GO Analysis
library(clusterProfiler)
library(stringr)
library(ggplot2)

# Convert gene symbols to Entrez IDs
gene_df <- bitr(diff_binded_genes, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
background_df <- bitr(background_genes, fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)
# remove duplicates
gene_entrez <- unique(gene_df$ENTREZID)
background_entrez <- unique(background_df$ENTREZID)

## Gene Ontology
go_bp <- enrichGO(gene = gene_entrez,
                  OrgDb = org.Hs.eg.db,
                  ont = "ALL",
                  keyType = "ENTREZID",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.2,
                  qvalueCutoff = 0.5,
                  readable = TRUE,
                  universe = background_entrez)

dotplot(go_bp, showCategory = 20) + ggtitle(NULL) + 
  theme(text=element_text(size=24), axis.text.x = element_text(size=24), axis.text.y = element_text(size=24), axis.title.x = element_text(size=24)) + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))




# ------------------
# Cross with RNAseq

res_NPC_KO_vs_WT <- read.csv('Results/DEGs/DEGs_NPC_KO_vs_WT.csv')
res_NPC_HET_vs_WT <- read.csv('Results/DEGs/DEGs_NPC_HET_vs_WT.csv')

NPC_KO_WT_sig_genes <- res_NPC_KO_vs_WT$X[res_NPC_KO_vs_WT$padj < 0.05]
NPC_HET_WT_sig_genes <- res_NPC_HET_vs_WT$X[res_NPC_HET_vs_WT$padj < 0.05]

DEFINE_COMPARISON <- 'KOvWT'
SEPARATE_REGULATION <-  'REPRESSED'           # 'REPRESSED'  or  'OPEN'   or   'NO'
FOLDCHANGE_THRESHOLD <- 0
PROMOTER_PEAKS <- T                # NOT implemented for the BOTH version

# Different code depending on if analyzing KOvWT, HETvWT, or combined KO+HETvWT.
# Also specify if separating by regulation


if(DEFINE_COMPARISON == 'BOTH'){
  print('BOTH KO+HET vs WT')
  res_deseq_WTvKO <- dba.report(db, method = DBA_DESEQ2, contrast = 2, th = 1)
  res_deseq_WTvHET <- dba.report(db, method = DBA_DESEQ2, contrast = 1, th = 1)
  
  annotation_matches_WTvKO <- findOverlaps(res_deseq_WTvKO, h38_gtf, maxgap = 5000)
  annotation_matches_WTvHET <- findOverlaps(res_deseq_WTvHET, h38_gtf, maxgap = 5000)
  res_deseq_WTvKO$SYMBOL <- NA_character_
  res_deseq_WTvKO$SYMBOL[queryHits(annotation_matches_WTvKO)] <- h38_gtf$gene_name[subjectHits(annotation_matches_WTvKO)]
  res_deseq_WTvHET$SYMBOL <- NA_character_
  res_deseq_WTvHET$SYMBOL[queryHits(annotation_matches_WTvHET)] <- h38_gtf$gene_name[subjectHits(annotation_matches_WTvHET)]
  
  diff_binded_genes_WTvKO <- unique(na.exclude(res_deseq_WTvKO$SYMBOL[res_deseq_WTvKO$FDR < 0.05]))
  diff_binded_genes_WTvHET <- unique(na.exclude(res_deseq_WTvHET$SYMBOL[res_deseq_WTvHET$FDR < 0.05]))
  diff_binded_genes_intersection <- intersect(diff_binded_genes_WTvKO, diff_binded_genes_WTvHET)
  
  intersected_genes <- intersect(diff_binded_genes_intersection, union(NPC_KO_WT_sig_genes,NPC_HET_WT_sig_genes))
  
  background_genes_intersection <- intersect(unique(na.exclude(res_deseq_WTvKO$SYMBOL)),unique(na.exclude(res_deseq_WTvHET$SYMBOL)))
  background_intersected_genes <- intersect(background_genes_intersection, res_NPC_KO_vs_WT$X)
  
} else if(DEFINE_COMPARISON == 'KOvWT'){
  if(SEPARATE_REGULATION == 'NO'){
    print('KOvWT, all intersections')
    intersected_genes <- intersect(diff_binded_genes, NPC_KO_WT_sig_genes)
    background_intersected_genes <- intersect(background_genes, res_NPC_KO_vs_WT$X)    # NOTE: diff_binded_genes needs updating depending on if looking at KO or HET
  } else if(SEPARATE_REGULATION == 'REPRESSED'){
    print('KOvWT, Repressed Genes')
    NPC_KO_WT_sig_genes <- res_NPC_KO_vs_WT$X[res_NPC_KO_vs_WT$padj < 0.05 & res_NPC_KO_vs_WT$log2FoldChange < -FOLDCHANGE_THRESHOLD]
    diff_binded_genes <- unique(na.exclude(res_deseq$SYMBOL[res_deseq$FDR < 0.05 & res_deseq$Fold < -FOLDCHANGE_THRESHOLD])) # Fold<0 means higher in Group 2 (ie KO)
    
    intersected_genes <- intersect(diff_binded_genes, NPC_KO_WT_sig_genes)
    background_intersected_genes <- intersect(background_genes, res_NPC_KO_vs_WT$X)
  } else if(SEPARATE_REGULATION == 'OPEN'){
    print('KOvWT, Open Genes')
    NPC_KO_WT_sig_genes <- res_NPC_KO_vs_WT$X[res_NPC_KO_vs_WT$padj < 0.05 & res_NPC_KO_vs_WT$log2FoldChange > FOLDCHANGE_THRESHOLD]
    diff_binded_genes <- unique(na.exclude(res_deseq$SYMBOL[res_deseq$FDR < 0.05 & res_deseq$Fold > FOLDCHANGE_THRESHOLD])) # Fold>0 means higher in Group 1 (ie WT)
    
    intersected_genes <- intersect(diff_binded_genes, NPC_KO_WT_sig_genes)
    background_intersected_genes <- intersect(background_genes, res_NPC_KO_vs_WT$X)
  }
} else if(DEFINE_COMPARISON == 'HETvWT'){
  if(SEPARATE_REGULATION == 'NO'){
    print('HETvWT, all intersections')
    intersected_genes <- intersect(diff_binded_genes, NPC_HET_WT_sig_genes)
    background_intersected_genes <- intersect(background_genes, res_NPC_HET_vs_WT$X)    # NOTE: diff_binded_genes needs updating depending on if looking at HET or HET
  } else if(SEPARATE_REGULATION == 'REPRESSED'){
    print('HETvWT, Repressed Genes')
    NPC_HET_WT_sig_genes <- res_NPC_HET_vs_WT$X[res_NPC_HET_vs_WT$padj < 0.05 & res_NPC_HET_vs_WT$log2FoldChange < -FOLDCHANGE_THRESHOLD]
    diff_binded_genes <- unique(na.exclude(res_deseq$SYMBOL[res_deseq$FDR < 0.05 & res_deseq$Fold < -FOLDCHANGE_THRESHOLD])) # Fold<0 means higher in Group 2 (ie HET)
    
    intersected_genes <- intersect(diff_binded_genes, NPC_HET_WT_sig_genes)
    background_intersected_genes <- intersect(background_genes, res_NPC_HET_vs_WT$X)
  } else if(SEPARATE_REGULATION == 'OPEN'){
    print('HETvWT, Open Genes')
    NPC_HET_WT_sig_genes <- res_NPC_HET_vs_WT$X[res_NPC_HET_vs_WT$padj < 0.05 & res_NPC_HET_vs_WT$log2FoldChange > FOLDCHANGE_THRESHOLD]
    diff_binded_genes <- unique(na.exclude(res_deseq$SYMBOL[res_deseq$FDR < 0.05 & res_deseq$Fold > FOLDCHANGE_THRESHOLD])) # Fold>0 means higher in Group 1 (ie WT)
    
    intersected_genes <- intersect(diff_binded_genes, NPC_HET_WT_sig_genes)
    background_intersected_genes <- intersect(background_genes, res_NPC_HET_vs_WT$X)
  }
}


# GO
gene_df <- bitr(intersected_genes, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
background_df <- bitr(background_intersected_genes, fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)

gene_entrez <- unique(gene_df$ENTREZID)
background_entrez <- unique(background_df$ENTREZID)

go_bp <- enrichGO(gene = gene_entrez,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  keyType = "ENTREZID",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.2,
                  qvalueCutoff = 0.4,
                  readable = TRUE,
                  universe = background_entrez)

dotplot(go_bp, showCategory = 20) + ggtitle(NULL) + 
  theme(text=element_text(size=24), axis.text.x = element_text(size=24), axis.text.y = element_text(size=24), axis.title.x = element_text(size=24)) + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))

# Extract genes related to term
GO_TERM_OF_INTEREST <- 'establishment of vesicle localization'

go_term_genes <- go_bp@result$geneID[go_bp@result$Description == GO_TERM_OF_INTEREST]

go_term_genes_parsed <- strsplit(go_term_genes, split = '/')[[1]]
writeClipboard(go_term_genes_parsed)












