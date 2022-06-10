# all required external package dependencies
# are installed with this script
# It should be the minimal set of packages
# required to run the code in 'lib.R' and 'run.R'
# ----- #
# option here to add an additional .libPaths()
# so that the "new" libraries are located in a
# isolated directory
# ---------------------- #
# Author: Stu Field
# ---------------------- #

# CRAN ----
install.packages("foreach")
install.packages("doParallel")
install.packages("tidyr")
install.packages("openxlsx")

install.packages("addGeneIDs")
install.packages("methAnalysis")
install.packages("getData")

# Bioconductor ----
install.packages("BiocManager")
cat("Using Bioconductor version:", as.character(BiocManager::version()), "\n")
BiocManager::install(version = "3.14")   # for R version 4.1
BiocManager::install('org.Bt.eg.db', force = TRUE)
BiocManager::install('clusterProfiler', force = TRUE)
BiocManager::install("TxDb.Btaurus.UCSC.bosTau9.refGene", force = TRUE)
BiocManager::install("BSgenome.Btaurus.UCSC.bosTau9", force = TRUE)
BiocManager::install("BSgenome", force = TRUE)
BiocManager::install(c("pathview", "gage", "gageData", "GenomicAlignments",
                       "TxDb.Hsapiens.UCSC.hg19.knownGene"), force = TRUE)
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", force = TRUE)
BiocManager::install("org.Hs.eg.db", force = TRUE)

BiocManager::install(c(
  "annotate", "AnnotationDbi", "backports", "ballgown", "BH", "Biobase",
  "BiocGenerics", "BiocParallel", "BiocVersion", "Biostrings", "bit",
  "blob", "DBI", "DelayedArray", "digest", "ellipsis", "genefilter",
  "GenomeInfoDb", "GenomeInfoDbData", "GenomicAlignments", "GenomicRanges",
  "IRanges", "lambda.r", "limma", "pkgconfig", "prettyunits", "Rcpp",
  "RCurl", "Rhtslib", "rlang", "Rsamtools", "RSQLite", "rtracklayer",
  "S4Vectors", "SummarizedExperiment", "sva", "vctrs", "XML", "XVector",
  "zlibbioc"), update = TRUE, ask = FALSE, force = TRUE)

BiocManager::install("methylKit", force = TRUE)
