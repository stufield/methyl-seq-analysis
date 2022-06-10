# This file contains all the functions required
# to run the analysis. It should be called from
# the 'run.R' file and should "set up" the
# desired analysis environment.
# No actual analysis is performed here, just the
# recipes for how to perform them.
# Objects and actual analysis should
# be generated in the 'run.R' file
# ------------------------------------- #
# Author: Stu Field
# ------------------------------------- #

library(genomation)
library(getData)

library(foreach)
library(doParallel)

library(data.table)
library(org.Bt.eg.db)
library(clusterProfiler)
library(BSgenome)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)
library(Biobase)
library(BSgenome.Btaurus.UCSC.bosTau9)
library(stringr)
library(tidyr)
library(pathview)
library(gage)
library(gageData)
library(GenomicAlignments)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(methylKit)


# analysis functions ----

# document func1
func1 <- function() {
  cat("great analysis!\n")
  invisible()
}

# document func2
func2 <- function() {
  cat("even better analysis!\n")
  invisible()
}
