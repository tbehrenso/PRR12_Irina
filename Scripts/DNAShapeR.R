library(DNAshapeR)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(fields)
library(slider)
library(tidyverse)
library(seqinr) # for read.fasta

MGW_THRESHOLD <- 4


ATAC_HEK <- read.table('/Users/tbehr/Desktop/HEK_ATAC/ATAC_HEK293.bed',header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="", skip = 1)
#ATAC_FASTA <- read.fasta('tmp.fa')
ATAC_FASTA <- read.fasta('/Users/tbehr/Desktop/HEK_ATAC/ATAC_HEK293_sequences.fa')

atac_granges <- GRanges(seqnames = ATAC_HEK$V1, ranges = IRanges(start = ATAC_HEK$V2, end = ATAC_HEK$V3))

# getFasta(atac_granges, BSgenome.Hsapiens.UCSC.hg38, width = 200, filename = "tmp.fa")


#pred <- getShape('tmp.fa')
pred <- getShape('/Users/tbehr/Desktop/HEK_ATAC/ATAC_HEK293_sequences.fa')

plotShape(pred$MGW)

# heatShape(pred$MGW, 20) --> gives error with custom fasta


# MGW smoothing function
mgw_windowed <- apply(pred$MGW, 1, function(x) slide_dbl(x, mean, .before = 4, .after = 4))


mgw_windowed_filtered <- mgw_windowed[,which(colSums(mgw_windowed < MGW_THRESHOLD, na.rm = T) > 1)]


# ----- Get AT content ---------
library(Biostrings)

# dna <- readDNAStringSet("tmp.fa")
dna <- readDNAStringSet('/Users/tbehr/Desktop/HEK_ATAC/ATAC_HEK293_sequences.fa')

# sliding AT content
win <- 9
win_wing <- ((win-1)/2)

at_content <- lapply(dna, function(x) {
  v <- letterFrequencyInSlidingView(x, view.width = win, letters = c("A","T"), as.prob = T)
  vsum <- rowSums(v)
  c(rep(NA, win_wing), vsum, rep(NA, win_wing))
})

at_content_filtered <- at_content[which(colSums(mgw_windowed < MGW_THRESHOLD, na.rm = T) > 1)]

at_content_filtered_smoothed <- lapply(at_content_filtered, function(x) slide_dbl(x, mean, .before = 4, .after = 4))

# ------ Filter by AT Content -------


# Convert the list of per-sequence AT-content vectors into a matrix
at_matrix <- do.call(cbind, at_content_filtered_smoothed)

stopifnot(identical(dim(at_matrix), dim(mgw_windowed_filtered)))

AT_THRESHOLD <- 0.65

# Elementwise: TRUE where a position passes BOTH filters simultaneously
combined_pass <- (at_matrix > AT_THRESHOLD) & (mgw_windowed_filtered < MGW_THRESHOLD)

# Which sequences have at least one position passing both?
passing_cols <- which(colSums(combined_pass, na.rm = TRUE) > 0)

# Map back to the actual sequences (since dna/ATAC_HEK were filtered the same way)
mgw_pass_idx <- which(colSums(mgw_windowed < MGW_THRESHOLD, na.rm = TRUE) > 1)
final_sequences <- dna[mgw_pass_idx][passing_cols]
final_bed_rows  <- ATAC_HEK[mgw_pass_idx, ][passing_cols, ]


# ------ Find consecutive TRUEs ------- (taken from Claude)

MIN_RUN <- 4  # tune based on expected footprint size

get_max_run <- function(col) {
  r <- rle(col)
  if (!any(r$values, na.rm = TRUE)) return(0)
  max(r$lengths[r$values], na.rm = TRUE)
}

run_lengths <- apply(combined_pass, 2, get_max_run)
strong_candidates <- which(run_lengths >= MIN_RUN)


strong_candidate_sequences <- dna[mgw_pass_idx][strong_candidates]
strong_candidate_bed_rows  <- ATAC_HEK[mgw_pass_idx, ][strong_candidates, ]

strong_candidate_grange <- GRanges(seqnames = strong_candidate_bed_rows$V1, ranges = IRanges(start = strong_candidate_bed_rows$V2, end = strong_candidate_bed_rows$V3))

# export.bed(strong_candidate_grange, 'strong_candidates.bed')


# ------ Cross with PRR12 peaks -------
process_narrowpeak <- function(peakfile_path){
  path_normalized <- normalizePath(peakfile_path)
  peakfile_read <- read.csv(peakfile_path, sep = '\t', header=F, col.names = c('chrom','chromStart','chromEnd','name','score','strand','signalValue','pValue','qValue','peak'))
  peakfile_read$strand <- '*'
  peak_range <- GRanges(seqnames=peakfile_read$chrom, ranges = IRanges(start = peakfile_read$chromStart, end = peakfile_read$chromEnd),
                        score = peakfile_read$signalValue, strand = peakfile_read$strand, signalValue = peakfile_read$signalValue, pValue = peakfile_read$pValue, qValue = peakfile_read$qValue, peak = peakfile_read$peak)
}

prr12_grange <- process_narrowpeak('Results/CUTandTag/External/GSM8590159_THC0467WNHN_PRR12_ChIP2_hg38_Hamed_peaks.narrowPeak')
prr12_pfilt8_grange <- prr12_grange[prr12_grange$qValue > 8,]


mgw_grange <- GRanges(seqnames = final_bed_rows$V1, ranges = IRanges(start = final_bed_rows$V2, end = final_bed_rows$V3))

prr12_mgw_overlap_hits <- findOverlaps(prr12_pfilt8_grange, mgw_grange, ignore.strand=T)

prr12_mgw_overlap <- pintersect(
  prr12_pfilt8_grange[queryHits(prr12_mgw_overlap_hits)],
  mgw_grange[subjectHits(prr12_mgw_overlap_hits)]
)

# export.bed(mgw_grange, 'mgw_at_sites.bed')



# ------ Filter by Canonical binding motif -------

library(Biostrings)
core_pattern <- "AATT"  # alternative to try: AAWWTT

has_core_motif <- sapply(dna[mgw_pass_idx], function(seq) {
  countPattern(core_pattern, seq, max.mismatch = 0) > 0
})

AATT_sequences <- dna[mgw_pass_idx][has_core_motif]
AATT_bed_rows  <- ATAC_HEK[mgw_pass_idx, ][has_core_motif, ]

AATT_grange <- GRanges(seqnames = AATT_bed_rows$V1, ranges = IRanges(start = AATT_bed_rows$V2, end = AATT_bed_rows$V3))

# export.bed(AATT_grange, 'AATT_motifs.bed')

# ------ Filter by Canonical binding motif, restricted to the passing window -------
library(Biostrings)
core_pattern <- "AATT"  # consider fixed = FALSE with an IUPAC pattern like "AAWWTT"
PAD <- 3                # small buffer around the run, since smoothing can shift edges

# For each strong candidate, find the start/end of its longest passing run
get_longest_run_range <- function(col) {
  r <- rle(col)
  vals <- r$values
  vals[is.na(vals)] <- FALSE  # treat NA as non-passing
  if (!any(vals)) return(c(NA_integer_, NA_integer_))
  true_idx <- which(vals)
  best <- true_idx[which.max(r$lengths[true_idx])]
  end_pos <- cumsum(r$lengths)[best]
  start_pos <- end_pos - r$lengths[best] + 1
  c(start_pos, end_pos)
}

run_ranges <- apply(combined_pass[, strong_candidates, drop = FALSE], 2, get_longest_run_range)
# run_ranges: row 1 = start, row 2 = end, one column per strong candidate

has_core_motif_in_window <- mapply(function(i, start, end) {
  if (is.na(start)) return(FALSE)
  seq_obj <- strong_candidate_sequences[[i]]
  win_start <- max(1, start - PAD)
  win_end   <- min(length(seq_obj), end + PAD)
  sub_seq   <- subseq(seq_obj, start = win_start, end = win_end)
  countPattern(core_pattern, sub_seq, max.mismatch = 0) > 0
}, seq_along(strong_candidate_sequences), run_ranges[1, ], run_ranges[2, ])

AATT_sequences <- strong_candidate_sequences[has_core_motif_in_window]
AATT_bed_rows  <- strong_candidate_bed_rows[has_core_motif_in_window, ]

AATT_grange <- GRanges(seqnames = AATT_bed_rows$V1,
                       ranges = IRanges(start = AATT_bed_rows$V2, end = AATT_bed_rows$V3))

# export.bed(AATT_grange, 'AATT_motifs.bed')





























