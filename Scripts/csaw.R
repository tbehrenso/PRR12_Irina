library(csaw)

COMPARISON <- 'WT-KO'


design_df <- data.frame(Sample = c('WT1','WT2','WT3','HET1','HET2','HET3','KO1','KO2','KO3'),
                        Type = as.factor(rep(c('WT','HET','KO'),each=3)))

design <- model.matrix(~design_df$Type)
rownames(design) <- design_df$Sample
colnames(design) <- c('WT','HET','KO')

# Load BAMs
bam.files <- c(
  '/beegfs/scratch/ric.broccoli/behrens.thomas/PRR12_Irina/CUTAndTag/H3k27ac/nf_out_macs2_broad_normalized_narrow/02_alignment/bowtie2/target/markdup/WT_R1.target.markdup.sorted.bam',
  '/beegfs/scratch/ric.broccoli/behrens.thomas/PRR12_Irina/CUTAndTag/H3k27ac/nf_out_macs2_broad_normalized_narrow/02_alignment/bowtie2/target/markdup/WT_R2.target.markdup.sorted.bam',
  '/beegfs/scratch/ric.broccoli/behrens.thomas/PRR12_Irina/CUTAndTag/H3k27ac/nf_out_macs2_broad_normalized_narrow/02_alignment/bowtie2/target/markdup/WT_R3.target.markdup.sorted.bam',
  '/beegfs/scratch/ric.broccoli/behrens.thomas/PRR12_Irina/CUTAndTag/H3k27ac/nf_out_macs2_broad_normalized_narrow/02_alignment/bowtie2/target/markdup/HET_R1.target.markdup.sorted.bam',
  '/beegfs/scratch/ric.broccoli/behrens.thomas/PRR12_Irina/CUTAndTag/H3k27ac/nf_out_macs2_broad_normalized_narrow/02_alignment/bowtie2/target/markdup/HET_R2.target.markdup.sorted.bam',
  '/beegfs/scratch/ric.broccoli/behrens.thomas/PRR12_Irina/CUTAndTag/H3k27ac/nf_out_macs2_broad_normalized_narrow/02_alignment/bowtie2/target/markdup/HET_R3.target.markdup.sorted.bam',
  '/beegfs/scratch/ric.broccoli/behrens.thomas/PRR12_Irina/CUTAndTag/H3k27ac/nf_out_macs2_broad_normalized_narrow/02_alignment/bowtie2/target/markdup/KO_R1.target.markdup.sorted.bam',
  '/beegfs/scratch/ric.broccoli/behrens.thomas/PRR12_Irina/CUTAndTag/H3k27ac/nf_out_macs2_broad_normalized_narrow/02_alignment/bowtie2/target/markdup/KO_R2.target.markdup.sorted.bam',
  '/beegfs/scratch/ric.broccoli/behrens.thomas/PRR12_Irina/CUTAndTag/H3k27ac/nf_out_macs2_broad_normalized_narrow/02_alignment/bowtie2/target/markdup/KO_R3.target.markdup.sorted.bam'
)

if(COMPARISON=='WT-KO'){
  design_df <- data.frame(Sample = c('WT1','WT2','WT3','KO1','KO2','KO3'),
                          Type = as.factor(rep(c('WT','KO'),each=3)))
  design <- model.matrix(~design_df$Type)
  rownames(design) <- design_df$Sample
  colnames(design) <- c('WT','KO')
  
  bam.files <- bam.files[-c(4,5,6)]
}else if(COMPARISON=='WT-HET'){
  design_df <- data.frame(Sample = c('WT1','WT2','WT3','HET1','HET2','HET3'),
                          Type = as.factor(rep(c('WT','HET'),each=3)))
  design <- model.matrix(~design_df$Type)
  rownames(design) <- design_df$Sample
  colnames(design) <- c('WT','HET')
  
  bam.files <- bam.files[-c(7,8,9)]
} else if(COMPARISON=='HET-KO'){
  design_df <- data.frame(Sample = c('HET1','HET2','HET3','KO1','KO2','KO3'),
                          Type = as.factor(rep(c('HET','KO'),each=3)))
  design <- model.matrix(~design_df$Type)
  rownames(design) <- design_df$Sample
  colnames(design) <- c('HET','KO')
  
  bam.files <- bam.files[-c(1,2,3)]
}


param <- readParam(minq=10, pe='both')
data <- windowCounts(bam.files, ext=100, param=param)


# Filter out uninteresting regions
binned <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)
keep <- filterWindowsGlobal(data, binned)$filter > log2(2)
summary(keep)
data <- data[keep,]

# Get the background for GO
avg.cpm.all <- rowMeans(cpm(asDGEList(data)))

keep.cpm <- avg.cpm.all > 0.5
data.universe <- data[keep.cpm, ]

universe_ranges <- rowRanges(data.universe)
saveRDS(universe_ranges, 'universe_ranges.RData')

# NEW -- Filter by local enrichment
surrounds <- 2000
neighbor <- suppressWarnings(resize(rowRanges(data), surrounds, fix="center"))
wider <- regionCounts(bam.files, regions=neighbor, ext=100, param=param)


# Calc normalization factors
data <- normFactors(binned, se.out=data)

# Identify differentially bound windows
library(edgeR)
y <- asDGEList(data)
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design, robust=TRUE)
results <- glmQLFTest(fit)

# Correct for multiple testing
merged <- mergeResults(data, results$table, tol=1000L, merge.args=list(max.width = 20000L))

clustered <- clusterWindows(rowRanges(data), rowData(data),
                            target=0.05, tol=1000L)


saveRDS(merged, 'merged.RData')




