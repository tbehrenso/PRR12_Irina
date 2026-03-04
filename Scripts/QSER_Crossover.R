library(tidyverse)
library(readxl)
library(GenomicRanges)

dixon_supp1 <- read_xlsx("/Users/tbehr/Desktop/NIHMS1701630-supplement-Data_S1.xlsx", sheet = 'QSER1 KO') %>%
  separate(`Region (hg38)`, into = c("chr", "coords"), sep = ":", remove=F) %>%
  separate(coords, into = c("start", "end"), sep = "-", convert = TRUE)


dixon_gr <- makeGRangesFromDataFrame(dixon_supp1, seqnames.field = 'chr', keep.extra.columns = T)

# (object from DiffBind.R)
res_deseq_sig <- res_deseq[res_deseq$FDR<0.05,]


# Intersect granges

grange_overlap <- findOverlaps(dixon_gr, res_deseq_sig, maxgap = 500, ignore.strand=T, select='first')

unique(res_deseq_sig$gene_id[na.exclude(grange_overlap)])
res_temp <- res_deseq_sig[na.exclude(grange_overlap),]

unique(res_temp$gene_id[res_temp$Fold > 0])
unique(res_temp$gene_id[res_temp$Fold < 0])
