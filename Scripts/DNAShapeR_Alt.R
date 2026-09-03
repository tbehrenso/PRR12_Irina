library(DNAshapeR)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(fields)
library(slider)
library(tidyverse)
library(seqinr)
library(Biostrings)
library(rtracklayer)

# ---- Parameters ----
MGW_THRESHOLD <- 4
AT_THRESHOLD  <- 0.65
MIN_RUN       <- 3
PAD           <- 3
win           <- 9
win_wing      <- (win - 1) / 2

# ---- Load data ----
ATAC_HEK <- read.table('/Users/tbehr/Desktop/HEK_ATAC/ATAC_HEK293.bed',
                       header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                       quote = "", skip = 1)

atac_granges <- GRanges(seqnames = ATAC_HEK$V1,
                        ranges   = IRanges(start = ATAC_HEK$V2 + 1, end = ATAC_HEK$V3))

# ---- DNAshapeR: always fixed 200bp windows ----
resized_granges <- resize(atac_granges, width = 200, fix = "center")

# getFasta(resized_granges, BSgenome.Hsapiens.UCSC.hg38, width = 200, filename = "tmp.fa")
pred <- getShape('tmp.fa', shapeType = 'MGW')   # reliably nSeq x 200

stopifnot(nrow(pred$MGW) == length(atac_granges))

plotShape(pred$MGW)
heatShape(pred$MGW, 20)

# smooth MGW across the 200bp window, one vector per sequence
mgw_windowed_list <- lapply(seq_len(nrow(pred$MGW)), function(i) {
  slide_dbl(pred$MGW[i, ], mean, .before = 4, .after = 4, na.rm=T)
})

# ---- AT content: real variable-length peak sequences ----
dna_peaks <- readDNAStringSet('/Users/tbehr/Desktop/HEK_ATAC/ATAC_HEK293_sequences.fa')

stopifnot(length(dna_peaks) == length(atac_granges))

at_content_smoothed <- lapply(dna_peaks, function(x) {
  v    <- letterFrequencyInSlidingView(x, view.width = win,
                                       letters = c("A", "T"), as.prob = TRUE)
  padded <- c(rep(NA, win_wing), rowSums(v), rep(NA, win_wing))
  padded
  #slide_dbl(padded, mean, .before = 4, .after = 4, na.rm=T)                         # remove this line (and added 'padded' as output above) to see if smoothing was losing info
})

# confirm AT vectors are peak-length (one value per bp)
stopifnot(all(lengths(at_content_smoothed) == width(atac_granges)))

# ---- Align MGW (window-space) and AT content (peak-space) ----
# Compute genomic overlap between each peak and its 200bp shape window.
# This correctly handles peaks both shorter AND longer than 200bp:
#   shorter peaks: window overhangs on both sides, overlap = peak itself
#   longer peaks:  peak overhangs on both sides, overlap = window itself
overlap_start_g <- pmax(start(atac_granges), start(resized_granges))
overlap_end_g   <- pmin(end(atac_granges),   end(resized_granges))
has_overlap     <- overlap_start_g <= overlap_end_g

# 1-based positions within the AT content vector (peak-space)
at_ov_start  <- overlap_start_g - start(atac_granges) + 1L
at_ov_end    <- overlap_end_g   - start(atac_granges) + 1L

# 1-based positions within the MGW vector (200bp-window-space)
mgw_ov_start <- overlap_start_g - start(resized_granges) + 1L
mgw_ov_end   <- overlap_end_g   - start(resized_granges) + 1L

# Build combined_pass in peak-coordinate space.
# Positions outside the overlap region (i.e. parts of long peaks beyond the
# 200bp window) are NA — they simply won't pass, which is correct since we
# have no shape information there.
combined_pass_list <- lapply(seq_along(at_content_smoothed), function(i) {
  result <- rep(NA, length(at_content_smoothed[[i]]))
  if (!has_overlap[i]) return(result)
  
  at_vec  <- at_content_smoothed[[i]]
  mgw_vec <- mgw_windowed_list[[i]]
  
  at_s  <- at_ov_start[i];  at_e  <- at_ov_end[i]
  mgw_s <- mgw_ov_start[i]; mgw_e <- mgw_ov_end[i]
  
  result[at_s:at_e] <- (at_vec[at_s:at_e] > AT_THRESHOLD) &
    (mgw_vec[mgw_s:mgw_e] < MGW_THRESHOLD)
  result
})

at_pass_list <- lapply(seq_along(at_content_smoothed), function(i) {
  result <- rep(NA, length(at_content_smoothed[[i]]))
  if (!has_overlap[i]) return(result)
  
  at_vec  <- at_content_smoothed[[i]]
  mgw_vec <- mgw_windowed_list[[i]]
  
  at_s  <- at_ov_start[i];  at_e  <- at_ov_end[i]
  mgw_s <- mgw_ov_start[i]; mgw_e <- mgw_ov_end[i]
  
  result[at_s:at_e] <- (at_vec[at_s:at_e] > AT_THRESHOLD)
  result
})

mgw_pass_list <- lapply(seq_along(at_content_smoothed), function(i) {
  result <- rep(NA, length(at_content_smoothed[[i]]))
  if (!has_overlap[i]) return(result)
  
  at_vec  <- at_content_smoothed[[i]]
  mgw_vec <- mgw_windowed_list[[i]]
  
  at_s  <- at_ov_start[i];  at_e  <- at_ov_end[i]
  mgw_s <- mgw_ov_start[i]; mgw_e <- mgw_ov_end[i]
  
  result[at_s:at_e] <- (mgw_vec[mgw_s:mgw_e] < MGW_THRESHOLD)
  result
})


# ---- Run-length filter ----
get_longest_run_range <- function(x) {
  x[is.na(x)] <- FALSE
  r <- rle(x)
  if (!any(r$values)) return(c(NA_integer_, NA_integer_))
  best     <- which(r$values)[which.max(r$lengths[which(r$values)])]
  end_pos  <- cumsum(r$lengths)[best]
  c(end_pos - r$lengths[best] + 1L, end_pos)
}

any_pass <- sapply(combined_pass_list, function(x) any(x, na.rm = TRUE))
at_pass <- sapply(at_pass_list, function(x) any(x, na.rm = TRUE))
mgw_pass <- sapply(mgw_pass_list, function(x) any(x, na.rm = TRUE))
pass_idx <- which(any_pass)

run_ranges_list <- lapply(combined_pass_list[pass_idx], get_longest_run_range)
run_lengths_vec <- sapply(run_ranges_list,
                          function(r) if (is.na(r[1])) 0L else r[2] - r[1] + 1L)

strong_mask       <- run_lengths_vec >= MIN_RUN
strong_global_idx <- pass_idx[strong_mask]   # indices back into atac_granges
run_ranges_strong <- run_ranges_list[strong_mask]

# ---- Motif search restricted to the passing window ----
# Searches only within the run (+ PAD) in peak-coordinate space,
# not across the full 200bp or full peak
core_pattern <- "AATT"  # use fixed = FALSE and e.g. "AAWWTT" for IUPAC degeneracy

has_core_motif <- mapply(function(gi, rr) {
  if (is.na(rr[1])) return(FALSE)
  seq_obj   <- dna_peaks[[gi]]
  win_start <- max(1L, rr[1] - PAD)
  win_end   <- min(length(seq_obj), rr[2] + PAD)
  sub_seq   <- subseq(seq_obj, start = win_start, end = win_end)
  countPattern(core_pattern, sub_seq, max.mismatch = 0, fixed = F) > 0
}, strong_global_idx, run_ranges_strong)

AATT_global_idx <- strong_global_idx[has_core_motif]
run_ranges_AATT <- run_ranges_strong[has_core_motif]

# ---- Build GRanges and export beds ----
# Helper: convert peak-space run ranges back to genomic coordinates.
# Works correctly for both short and long peaks since run_ranges are
# already expressed as offsets from the peak start (atac_granges).
make_subregion_granges <- function(global_idx, run_ranges) {
  base <- atac_granges[global_idx]
  GRanges(
    seqnames = seqnames(base),
    ranges   = IRanges(
      start = start(base) + sapply(run_ranges, `[`, 1) - 1L,
      end   = start(base) + sapply(run_ranges, `[`, 2) - 1L
    )
  )
}

# Tier 1: any position passing AT+MGW (least stringent)
passing_peaks_grange    <- atac_granges[pass_idx]

# Tier 2: passing run >= MIN_RUN (strong structural candidates)
strong_subregion_grange <- make_subregion_granges(strong_global_idx, run_ranges_strong)

# Tier 3: strong run + AATT motif within that run (highest confidence)
AATT_subregion_grange   <- make_subregion_granges(AATT_global_idx, run_ranges_AATT)

# export.bed(passing_peaks_grange,    'AT_MGW_passing_peaks.bed')
# export.bed(strong_subregion_grange, 'AT_MGW_strong_subregions.bed')
# export.bed(AATT_subregion_grange,   'AATT_subregions.bed')

# ---- PRR12 peak overlap ----
process_narrowpeak <- function(peakfile_path) {
  df <- read.csv(peakfile_path, sep = '\t', header = FALSE,
                 col.names = c('chrom','chromStart','chromEnd','name','score',
                               'strand','signalValue','pValue','qValue','peak'))
  df$strand <- '*'
  GRanges(seqnames    = df$chrom,
          ranges      = IRanges(start = df$chromStart, end = df$chromEnd),
          strand      = df$strand,
          score       = df$signalValue,
          signalValue = df$signalValue,
          pValue      = df$pValue,
          qValue      = df$qValue,
          peak        = df$peak)
}

filtered_grange_to_overlap <- passing_peaks_grange


prr12_grange        <- process_narrowpeak('Results/CUTandTag/External/GSM8590159_THC0467WNHN_PRR12_ChIP2_hg38_Hamed_peaks.narrowPeak')
prr12_pfilt8_grange <- prr12_grange[prr12_grange$qValue > 8]

# Overlap against the most stringent set (AATT-confirmed subregions)
prr12_pfilt8_overlap_hits <- findOverlaps(prr12_pfilt8_grange, filtered_grange_to_overlap,
                                   ignore.strand = TRUE)

prr12_pfilt8_overlap <- pintersect(
  prr12_pfilt8_grange[queryHits(prr12_pfilt8_overlap_hits)],
  filtered_grange_to_overlap[subjectHits(prr12_pfilt8_overlap_hits)],
  ignore.strand = TRUE
)

prr12_overlap_hits <- findOverlaps(prr12_grange, filtered_grange_to_overlap,
                                   ignore.strand = TRUE)
prr12_overlap <- pintersect(
  prr12_grange[queryHits(prr12_overlap_hits)],
  filtered_grange_to_overlap[subjectHits(prr12_overlap_hits)],
  ignore.strand = TRUE
)

# export.bed(prr12_overlap, 'PRR12_AATT_overlap.bed')

# ---- Annotate final ranges ----
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

prr12_pfilt8_overlap_peakAnno <- annotatePeak(
  prr12_pfilt8_overlap,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db",
  tssRegion = c(-3000, 3000)
)
prr12_overlap_peakAnno <- annotatePeak(
  prr12_overlap,
  TxDb = txdb,
  annoDb = "org.Hs.eg.db",
  tssRegion = c(-3000, 3000)
)

plotAnnoPie(prr12_overlap_peakAnno)

prr12_pfilt8_overlap_genes <- unique(prr12_pfilt8_overlap_peakAnno@anno$SYMBOL)
prr12_overlap_genes <- unique(prr12_overlap_peakAnno@anno$SYMBOL)

# Overlap with RNAseq
res_NPC_KO_vs_WT <- read.csv('Results/DEGs/DEGs_NPC_KO_vs_WT.csv')
res_NPC_HET_vs_WT <- read.csv('Results/DEGs/DEGs_NPC_HET_vs_WT.csv')

NPC_KO_WT_sig_genes_up <- na.exclude(res_NPC_KO_vs_WT$X[res_NPC_KO_vs_WT$padj < 0.05 & res_NPC_KO_vs_WT$log2FoldChange>1])  # UP in KO
NPC_KO_WT_sig_genes_down <- na.exclude(res_NPC_KO_vs_WT$X[res_NPC_KO_vs_WT$padj < 0.05 & res_NPC_KO_vs_WT$log2FoldChange< -1]) # DOWN in KO

intersect(prr12_overlap_genes, NPC_KO_WT_sig_genes_up)
intersect(prr12_overlap_genes, NPC_KO_WT_sig_genes_down)




at_pooled <- unlist(lapply(seq_along(at_content_smoothed), function(i) {
  if (!has_overlap[i]) return(NULL)
  at_content_smoothed[[i]][at_ov_start[i]:at_ov_end[i]]
}))

mgw_pooled <- unlist(lapply(seq_along(mgw_windowed_list), function(i) {
  if (!has_overlap[i]) return(NULL)
  mgw_windowed_list[[i]][mgw_ov_start[i]:mgw_ov_end[i]]
}))

hist(at_pooled)
hist(mgw_pooled)
cor(at_pooled, mgw_pooled, use = "pairwise.complete.obs")

set.seed(1)
shuffled_order <- sample(seq_along(mgw_windowed_list))

paired <- lapply(seq_along(at_content_smoothed), function(k) {
  i <- k                      # AT comes from peak i
  j <- shuffled_order[k]      # MGW comes from a randomly reassigned peak j
  if (!has_overlap[i] || !has_overlap[j]) return(NULL)
  at_vec  <- at_content_smoothed[[i]][at_ov_start[i]:at_ov_end[i]]
  mgw_vec <- mgw_windowed_list[[j]][mgw_ov_start[j]:mgw_ov_end[j]]
  n <- min(length(at_vec), length(mgw_vec))
  cbind(at = at_vec[1:n], mgw = mgw_vec[1:n])
})

paired_mat <- do.call(rbind, paired)
cor(paired_mat[,"at"], paired_mat[,"mgw"], use = "pairwise.complete.obs")



summary((mgw_windowed_list))
threshold_points <- c(3.5,3.75, 4, 4.25, 4.5, 4.75, 5)
threshold_points <- seq(3, 5, length.out = 15)
threshold_ratios <- sapply(threshold_points, function(t) mean(sapply(mgw_windowed_list, function(v) any(v < t, na.rm=TRUE))))
threshold_ratios

plot(x=threshold_points, y=threshold_ratios, type='l')






