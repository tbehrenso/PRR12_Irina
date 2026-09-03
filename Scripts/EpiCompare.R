library(EpiCompare)
library(GenomicRanges)


process_narrowpeak <- function(peakfile_path){
  path_normalized <- normalizePath(peakfile_path)
  peakfile_read <- read.csv(peakfile_path, sep = '\t', header=F, col.names = c('chrom','chromStart','chromEnd','name','score','strand','signalValue','pValue','qValue','peak'))
  peakfile_read$strand <- '*'
  peak_range <- GRanges(seqnames=peakfile_read$chrom, ranges = IRanges(start = peakfile_read$chromStart, end = peakfile_read$chromEnd),
          score = peakfile_read$signalValue, strand = peakfile_read$strand, signalValue = peakfile_read$signalValue, pValue = peakfile_read$pValue, qValue = peakfile_read$qValue, peak = peakfile_read$peak)
}

prr12_grange <- process_narrowpeak('Results/CUTandTag/External/GSM8590159_THC0467WNHN_PRR12_ChIP2_hg38_Hamed_peaks.narrowPeak')
prr12_pfilt8_grange <- prr12_grange[prr12_grange$qValue > 8,]
h3k4me1_grange <- process_narrowpeak('Results/CUTandTag/External/ENCODE_or_ChIPAtlas/h3k4me1_ENCFF351MDO.bed')
h3k27ac_grange <- process_narrowpeak('Results/CUTandTag/External/ENCODE_or_ChIPAtlas/h3k27ac_ENCFF225FPN.bed')
ctcf_grange <- process_narrowpeak('Results/CUTandTag/External/ENCODE_or_ChIPAtlas/CTCF_ENCFF521EHU.bed')
h3k4me3_grange <- process_narrowpeak('Results/CUTandTag/External/ENCODE_or_ChIPAtlas/h3k4me3_ENCFF617NUV.bed')
h3k9ac_grange <- process_narrowpeak('Results/CUTandTag/External/ENCODE_or_ChIPAtlas/h3k9ac_ENCFF138PYK_misc.bed')
h3k9me3_grange <- process_narrowpeak('Results/CUTandTag/External/ENCODE_or_ChIPAtlas/h3k9me3_ENCFF899MET.bed')
h3k36me3_grange <- process_narrowpeak('Results/CUTandTag/External/ENCODE_or_ChIPAtlas/h3k36me3_ENCFF655QDU.bed')
h3k4me1_chip_grange <- process_narrowpeak('Results/CUTandTag/External/ENCODE_or_ChIPAtlas/h3k4me1_ENCFF301UTR.bed')
setdb1_grange <- process_narrowpeak('Results/CUTandTag/External/ENCODE_or_ChIPAtlas/SETDB1_ENCFF812YDP.bed')
trim28_grange <- process_narrowpeak('Results/CUTandTag/External/ENCODE_or_ChIPAtlas/TRIM28_ENCFF032VPU.bed')
znf263_grange <- process_narrowpeak('Results/CUTandTag/External/ENCODE_or_ChIPAtlas/ZNF263_ENCFF108NXV.bed')
h3k27me3_bed12 <- read.csv('Results/CUTandTag/External/ENCODE_or_ChIPAtlas/h3k27me3_external_c1.5_C1.00_l200_g30_G800_broad.bed12', sep = '\t', header=F,
         col.names = c('chrom','chromStart','chromEnd','name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts','signalValue','pValue','qValue'),skip = 1)
h3k27me3_bed12$strand <- '*'
h3k27me3_grange <- GRanges(seqnames=h3k27me3_bed12$chrom, ranges = IRanges(start = h3k27me3_bed12$chromStart, end = h3k27me3_bed12$chromEnd),
                           score = h3k27me3_bed12$score, strand = h3k27me3_bed12$strand, signalValue = h3k27me3_bed12$signalValue, pValue = h3k27me3_bed12$pValue, qValue = h3k27me3_bed12$qValue)
hela_h3k27me3_grange <- process_narrowpeak('Results/CUTandTag/External/ENCODE_or_ChIPAtlas/hela_h3k27me3_ENCFF252BLX.bed')
h3k27me3_custom_grange <- process_narrowpeak('Results/CUTandTag/External/ENCODE_or_ChIPAtlas/H3k27me3_GSE133391_peaks.narrowPeak')
rad21_grange <- process_narrowpeak('Results/CUTandTag/External/ENCODE_or_ChIPAtlas/RAD21_DRX013191.05.bed')


peakranges <- list(PRR12_ChIP = prr12_grange,
                   PRR12_ChIP_pfilt8 = prr12_pfilt8_grange,
                   H3K4ME1 = h3k4me1_grange,
                   H3K27AC = h3k27ac_grange,
                   CTCF = ctcf_grange,
                   H3K4ME3 = h3k4me3_grange,
                   H3k9AC = h3k9ac_grange,
                   H3K9ME3 = h3k9me3_grange,
                   H3k36ME3 = h3k36me3_grange)

peakranges <- list(PRR12_ChIP = prr12_grange,
                   PRR12_ChIP_pfilt8 = prr12_pfilt8_grange,
                   RAD21 = rad21_grange,
                   H3K27AC = h3k27ac_grange,
                   H3K4ME3 = h3k4me3_grange,
                   H3K9ME3 = h3k9me3_grange,
                   CTCF = ctcf_grange,
                   H3K27ME3_CUSTOM = h3k27me3_custom_grange
                   )

common_chroms <- Reduce(
  intersect,
  lapply(peakranges, function(x){
    unique(as.character(seqnames(x)))
  })
)

peaklist <- lapply(peakranges, function(x){
  x[as.character(seqnames(x)) %in% common_chroms]
})

# EpiCompare(peakfiles = peakranges,
#            genome_build = 'hg38',
#            genome_build_output = 'hg38',
#            output_dir = 'Results/CUTandTag/EpiCompare',
#            save_output=TRUE,
#            run_all = F, upset_plot = T, stat_plot = T, chromHMM_plot = T, chipseeker_plot = T, enrichment_plot = T, tss_plot = F, precision_recall_plot = T, corr_plot = T)

EpiCompare(peakfiles = peakranges[3:length(peakranges)],
           genome_build = 'hg38',
           genome_build_output = 'hg38',
           output_dir = 'Results/CUTandTag/EpiCompare',
           save_output=TRUE,
           reference = peakranges[2],
           chromHMM_annotation = 'H1hesc',
           run_all = F, upset_plot = T, stat_plot = T, chromHMM_plot = T, chipseeker_plot = T, enrichment_plot = T, tss_plot = F, precision_recall_plot = F, corr_plot = T)






# Extra: Analysis Overlap -------------------------------------------------

prr12_ctcf_overlap_hits <- findOverlaps(prr12_pfilt8_grange, ctcf_grange, ignore.strand=T)

prr12_ctcf_overlap <- pintersect(
  prr12_pfilt8_grange[queryHits(prr12_ctcf_overlap_hits)],
  ctcf_grange[subjectHits(prr12_ctcf_overlap_hits)]
)

plot_chromHMM(peaklist = list(prr12=prr12_pfilt8_grange, ctcf=ctcf_grange, overlap = prr12_ctcf_overlap),
              cell_line = 'H1hesc', 
              genome_build = 'hg38')















