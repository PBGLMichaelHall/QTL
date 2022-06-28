#!/usr/bin/env Rscript
# supress warnings
#options(warn=-1)


library(modeest)
library(ggplot2)
library(gtools)
library(dplyr)
library(readr)
library(tidyr)
library(Rcpp)
library(locfit)
library(knitr)
library(rmarkdown)
library(kableExtra)
library(devtools)
library(data.table)
library(vcfR)
#install.packages("optparse")
library(optparse)
devtools::install_github("pbgl/QTLseqr")
library(QTLseqr)


# option parsing #

parser <- OptionParser()

parser <- add_option(parser,
		     opt_str = c("-i", "--input"),
                     type = "character",
                     dest = 'vcf',
                     help = "Variant Called Format file or directory (required). Either a full path to a BSA vcf file or one of such files."
                     )

parser <- add_option(parser,
                     opt_str = c("-H", "--highBulk"),
                     type = "character",
                     dest = 'H',
                     help="Name of High Bulk Sample from BSA experiment"
                     )

parser <- add_option(parser,
                     opt_str = c("-l", "--lowBulk"),
                     type = "character",
                     dest = 'L',
                     help = "Name of Low Bulk Sample from BSA experiment"
                     )

parser <- add_option(parser,
                     opt_str = c("-c", "--chromosomes"),
                     type = "character",
                     dest = 'Chrom',
                     help = "Please Provide a text file with chromosome list in an R vector"
                     )

parser <- add_option(parser,
                     opt_str = c("-r", "--refAlleleFreq"),
                     type = "numeric",
                     dest = 'ref',
                     help = "Filtering SNP value for Reference Allele Frequency 0 < RAF < 1"
                     )

parser <- add_option(parser,
                     opt_str = c("-m", "--minTotalDepth"),
                     type = "numeric",
                     dest = 'minTD',
                     help = "Filtering SNP value for minimum Total Depth"
                     )

parser <- add_option(parser,
                     opt_str = c("-d", "--maxTotalDepth"),
                     type = "numeric",
                     dest = 'maxD',
                     help = "Filtering SNP value for maximum Total Depth"
                     )

parser <- add_option(parser,
                     opt_str = c("-D", "--minSampleDepth"),
                     type = "numeric",
                     dest = 'minSD',
                     help = "Filtering SNP value for minimum Sample Depth"
                     )

parser <- add_option(parser,
                     opt_str = c("-G", "--minGQ"),
                     type = "numeric",
                     dest = 'GQ',
                     help = "Filtering SNP value for minimum Genotype Quality"
                     )

parser <- add_option(parser,
                     opt_str = c("-w", "--windowSize"),
                     type = "numeric",
                     dest = "WS",
                     help = "Specify window size for counting SNPs etc."
                     )

parser <- add_option(parser,
                     opt_str = c("-f", "--filterThreshold"),
                     type = "numeric",
                     dest = "FT",
                     help = "Filter Threshold for more filtering of outliers"
                     )

parser <- add_option(parser,
                     opt_str = c("-P", "--popStruc"),
                     type = "character",
                     dest = "PS",
                     help = "Population Structure of BSA"
                     )

parser <- add_option(parser,
                     opt_str = c("-B", "--bulkSize"),
                     type = 'numeric',
                     dest = "BS",
                     help = "Provide a vector of bulk sizes c(Highbulk,LowBulk) ---. c(45,38)"
                     )

parser <- add_option(parser,
                     opt_str = c("-q","--qVAR"),
                     type = "numeric",
                     dest = "QQ",
                     help = "Provide a value for plotQTLStats q"
                     )

parser <- add_option(parser,
                     opt_str = c("-a","--alpha"),
                     type = "numeric",
                     dest = "alpha",
                     help = "Provide alpha level for QTLTABLe"
                     )

opt = parse_args(parser)




vcf = opt$vcf
H = opt$H
L = opt$L
Chrom = opt$Chrom
Chrom=scan("Chr.txt", character(), sep=",")

print(Chrom)

ref = opt$ref
minTD = opt$minTD
maxD = opt$maxD
minSD = opt$minSD
GQ = opt$GQ
WS = opt$WS
FT = opt$FT
PS = opt$PS
BS = opt$BS
BS=scan("bulksize.txt", character(), sep=",")
QQ = opt$QQ
AA = opt$alpha


df <-
    QTLseqr::importFromVCF(
        file = vcf,
        highBulk = H,
        lowBulk = L,
        chromList = Chrom
     )

#print(head(df,10))

df_filt <-
    filterSNPs(
        SNPset = df,
        refAlleleFreq = ref,
        minTotalDepth = minTD,
        maxTotalDepth = maxD,
        minSampleDepth = minSD,
        #depthDifference = 100,
        minGQ = GQ,
        verbose = TRUE
   )


#print(head(df_filt,10))

df_filt<-QTLseqr::runGprimeAnalysis(
  SNPset = df_filt,
  windowSize = WS,
  outlierFilter = "deltaSNP",
  filterThreshold = FT)

#print(head(df_filt,10))


df_filt <- QTLseqr::runQTLseqAnalysis(
  SNPset = df_filt,
  windowSize = WS,
  popStruc = PS,
  bulkSize = BS,
  replications = 10000,
  intervals = c(95, 99)
)

#print(head(df_filt))


pdf(file = "plotcentral.pdf")
QTLseqr::plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = QQ)
QTLseqr::plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals  = TRUE)
ggplot2::ggplot(data =df) + ggplot2::geom_histogram(ggplot2::aes(x = DP.LOW + DP.HIGH)) + ggplot2::xlim(0,400)
ggplot2::ggplot(data = df) + ggplot2::geom_histogram(ggplot2::aes(x = REF_FRQ))
ggplot2::ggplot(data = df) + ggplot2::geom_histogram(ggplot2::aes(x = DP.LOW))
ggplot2::ggplot(data = df) + ggplot2::geom_histogram(ggplot2::aes(x = DP.HIGH))
ggplot2::ggplot(data = df) + ggplot2::geom_histogram(ggplot2::aes(x = GQ.LOW))
ggplot2::ggplot(data = df) + ggplot2::geom_histogram(ggplot2::aes(x = GQ.HIGH))

dev.off()
#export summary CSV
QTLseqr::getQTLTable(SNPset = df_filt, alpha = AA, export = TRUE, fileName = "my_BSA_QTL.csv")
