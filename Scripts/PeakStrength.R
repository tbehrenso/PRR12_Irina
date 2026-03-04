


strength_matrix <- read.csv('/Users/tbehr/Desktop/CpG_H3K27me3_raw_unnormalized.tsv', sep='\t')
strength_matrix <- read.csv('/Users/tbehr/Desktop/CpG_H3K27me3_raw_unnormalized.fixed.tsv', sep='\t')
#strength_matrix <- read.csv('/Users/tbehr/Desktop/macs2_min3_H3K27me3_raw_unnormalized.tsv', sep='\t')


strength_matrix_noNA <- na.exclude(strength_matrix)


strength_matrix_noNA_reduced <- strength_matrix_noNA[,4:10]

colnames(strength_matrix_noNA_reduced) <- c('WT2','WT3','HET1','HET2','HET3','KO1','KO3')


coldata <- data.frame(
  condition = c('WT','WT','HET','HET','HET','KO','KO')
)
rownames(coldata) <- colnames(strength_matrix_noNA_reduced)

# Normalize??
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = round(strength_matrix_noNA_reduced),
  colData = coldata,
  design = ~ condition
)

dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized=TRUE)


# strength_matrix_noNA_reduced <- norm_counts

# Compute Peak Strength
WT_mean <- rowMeans(strength_matrix_noNA_reduced[, coldata$condition == 'WT'])
HET_mean <- rowMeans(strength_matrix_noNA_reduced[, coldata$condition == 'HET'])
KO_mean <- rowMeans(strength_matrix_noNA_reduced[, coldata$condition == 'KO'])
WT_logmean <- log(WT_mean+1)
HET_logmean <- log(HET_mean+1)
KO_logmean <- log(KO_mean+1)


boxplot(
  WT_mean, HET_mean, KO_mean,
  names = c("WT", 'HET',"KO"),
  ylab = "H3K27me3 signal",
  main = "H3K27me3 at CpG islands",
  ylim = c(0,600)
)
boxplot(
  WT_logmean, HET_logmean, KO_logmean,
  names = c("WT", 'HET',"KO"),
  ylab = "H3K27me3 signal",
  main = "H3K27me3 at CpG islands"
)

# logFC
log2FC_KO_WT <- log2((KO_mean + 1) / (WT_mean + 1))
log2FC_HET_WT <- log2((HET_mean + 1) / (WT_mean + 1))
log2FC_KO_HET <- log2((KO_mean + 1) / (HET_mean + 1))

hist(
  log2FC_KO_WT,
  breaks = 100,
  main = "Log2(KO / WT) H3K27me3 at CpG islands",
  xlab = "log2 fold change"
)
hist(
  log2FC_HET_WT,
  breaks = 100,
  main = "Log2(HET / WT) H3K27me3 at CpG islands",
  xlab = "log2 fold change"
)
hist(
  log2FC_KO_HET,
  breaks = 100,
  main = "Log2(KO / HET) H3K27me3 at CpG islands",
  xlab = "log2 fold change"
)

mean(KO_mean)
mean(HET_mean)
mean(WT_mean)

t.test(KO_logmean, y=WT_logmean, paired=T)
t.test(HET_logmean, y=WT_logmean, paired=T)
t.test(KO_logmean, y=HET_logmean, paired=T)

# Cohen's d
mean(KO_mean - WT_mean) / sd(KO_mean - WT_mean)
mean(HET_mean - WT_mean) / sd(HET_mean - WT_mean)
mean(KO_mean - HET_mean) / sd(KO_mean - HET_mean)



# ---------------------------------------------------------------
# Scaling Factor
# based on https://www.biostars.org/p/413626/#414440
# ---------------------------------------------------------------
library(edgeR)

strength_matrix_noNA_reduced

## edgeR:: calcNormFactors
NormFactor <- calcNormFactors(object = strength_matrix_noNA_reduced, method = "TMM")
# raw library size
LibSize <- colSums(strength_matrix_noNA_reduced)
# calc size factors
SizeFactors <- NormFactor * LibSize / 1000000

# reciprocal size factors
SizeFactors.Reciprocal <- 1/SizeFactors


# ---------------------------------------------------------------
# Bar plots
# ---------------------------------------------------------------
library(ggplot2)

samplenames <- c('WT1','WT2','WT3','HET1','HET2','HET3','KO1','KO2','KO3')
bp_sum_macs2_broad <- c(12902056,5431347,7826837,4469304,2022320,1348523,1841080,2500449,1455280)
peak_count_macs2_broad <- c(13629,8282,15030,6549,4614,2891,4573,9631,3351)

samplenames <- c('WT2','WT3','HET1','HET2','HET3','KO1','KO3')
sampletype <- c('WT','WT','HET','HET','HET','KO','KO')
bp_sum_macs2_broad <- c(5431347,7826837,4469304,2022320,1348523,1841080,1455280)
peak_count_macs2_broad <- c(8282,15030,6549,4614,2891,4573,3351)



barplot_df <- data.frame(
  samplename = samplenames,
  sampletype = factor(sampletype, levels = c('WT','HET','KO')),
  bp_sum_macs2_broad = bp_sum_macs2_broad,
  peak_count_macs2_broad = peak_count_macs2_broad
)

ggplot(barplot_df, aes(x = sampletype, y = peak_count_macs2_broad)) +  # can change y depending on focus
  stat_summary(
    fun = mean,
    geom = "col",
    fill = "grey80",
    color = "black",
    width = 0.6
  ) +
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  theme_classic() +
  labs(
    x = "Condition",
    y = "Peak Count" 
  )



# scatterplot with y=x line

condition_mean_df <- data.frame(WT = WT_mean, KO = KO_mean, HET = HET_mean)

ggplot(condition_mean_df, aes(x = WT, y = KO)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype='dashed', col='red') +
  theme_minimal() +
  labs(x='WT',y='KO') +
  coord_equal()







