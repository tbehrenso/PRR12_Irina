library(LOLA)
library(GenomicRanges)

process_narrowpeak <- function(peakfile_path){
  path_normalized <- normalizePath(peakfile_path)
  peakfile_read <- read.csv(peakfile_path, sep = '\t', header=F, col.names = c('chrom','chromStart','chromEnd','name','score','strand','signalValue','pValue','qValue','peak'))
  peakfile_read$strand <- '*'
  peak_range <- GRanges(seqnames=peakfile_read$chrom, ranges = IRanges(start = peakfile_read$chromStart, end = peakfile_read$chromEnd),
                        score = peakfile_read$signalValue, strand = peakfile_read$strand, signalValue = peakfile_read$signalValue, pValue = peakfile_read$pValue, qValue = peakfile_read$qValue, peak = peakfile_read$peak)
}

# Calculate TAD Boundaries ------------------------------------------------
tad_bed <- readBed('data/External/TADMap_scaffold_hs.bed')

tads <- sort(tad_bed)

boundary.pos <- end(tads[-length(tads)])

boundaries <- GRanges(
  seqnames = seqnames(tads[-length(tads)]),
  ranges = IRanges(boundary.pos, width = 1)
)

boundaries10kb <- resize(
  boundaries,
  width = 10001,
  fix = "center"
)

# LOLA --------------------------------------------------------------------


# Load Files
regionDB = loadRegionDB("LOLACore/hg38")

tad_bed <- readBed('data/External/TADMap_scaffold_hs.bed')
prr12_bed <- readBed('Results/CUTandTag/External/GSM8590159_THC0467WNHN_PRR12_ChIP2_hg38_Hamed_peaks.narrowPeak')

prr12_grange <- process_narrowpeak('Results/CUTandTag/External/GSM8590159_THC0467WNHN_PRR12_ChIP2_hg38_Hamed_peaks.narrowPeak')
prr12_pfilt8_grange <- prr12_grange[prr12_grange$qValue > 8,]

# LOLA RDS Files
hg38_permissive <- readRDS('data/External/LOLA/hg38_permissive_UniBind_LOLA.RDS')
hg38_permissive_universe <- readRDS('data/External/LOLA/hg38_permissive_UniBind_LOLA_universe.RDS')
hg38_robust <- readRDS('data/External/LOLA/hg38_robust_UniBind_LOLA.RDS')
hg38_robust_universe <- readRDS('data/External/LOLA/hg38_robust_UniBind_LOLA_universe.RDS')

userSets = GRangesList(prr12_pfilt8_grange)
locResults = runLOLA(userSets, boundaries10kb, hg38_permissive, cores=1)



# regioneR ----------------------------------------------------------------
library(regioneR)
library(BSgenome.Hsapiens.UCSC.hg38)

genome_hg38 <- BSgenome.Hsapiens.UCSC.hg38

numOverlaps(boundaries10kb, prr12_pfilt8_grange, count.once=TRUE)

pt <- overlapPermTest(A=boundaries10kb, B=prr12_pfilt8_grange, ntimes=50, genome=genome_hg38)

plot(pt)















