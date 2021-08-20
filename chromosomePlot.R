# ===== Install Stuff =====

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install("karyoploteR")

BiocManager::install("BSgenome", force=TRUE)

library(BiocManager)
install("BSgenome.Hsapiens.UCSC.hg38")

library(karyoploteR)

# ===== Generate Plot =====

regions1 <- regioneR::toGRanges("chromDataForPlot/miscalls_grch38_bwameth_emseq_200ng_4cycles_LB_combined_reps_bed.sorted.bed")

regions2 <- regioneR::toGRanges("chromDataForPlot/grch38_200ng_bwameth_emseq_cleaned_cut.tsv")

chromsToPlot <- c("chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

# chromsToPlot <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12")

chromDataX1 = list()

chromDataY1 = list()

chromDataX2 = list()

chromDataY2 = list()

i <- 1

for(chrom in chromsToPlot)
{
  
kp.tmp <- karyoploteR::plotKaryotype("hg38", chr=chrom)

kp.tmp <- karyoploteR::kpPlotDensity(kp.tmp, data=regions1)

density1 <- kp.tmp$latest.plot$computed.values$density
windows1 <- kp.tmp$latest.plot$computed.values$windows

density1[density1>30000] <- 30000

kp.tmp <- karyoploteR::kpPlotDensity(kp.tmp, data=regions2)

density2 <- kp.tmp$latest.plot$computed.values$density
windows2 <- kp.tmp$latest.plot$computed.values$windows

density2[density2>30000] <- 30000

chromDataY1[[i]] = density1
chromDataX1[[i]] = windows1

chromDataY2[[i]] = density2
chromDataX2[[i]] = windows2

i <- i+1

}

i <- 1

kp <- karyoploteR::plotKaryotype("hg38", chr=chromsToPlot)

for(chrom in chromsToPlot)
{
  karyoploteR::kpArea(kp, chromDataX2[[i]], y=chromDataY2[[i]], chr=chrom, col="blue", ymax=max(30000), ymin=0, r0=0.0, r1=0.5)
  karyoploteR::kpArea(kp, chromDataX1[[i]], y=chromDataY1[[i]], chr=chrom, col="red", ymax=max(30000), ymin=0, r0=0.5, r1=1.0)
  i <- i+1
}