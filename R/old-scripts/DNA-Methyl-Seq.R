install.packages("foreach")
install.packages("doParallel")
install.packages("BiocManager")
install.packages("tidyr")
install.packages("openxlsx")
BiocManager::install("BiocGenerics", force = TRUE)

install.packages("addGeneIDs")
install.packages("methAnalysis")
install.packages("getData")

BiocManager::install('org.Bt.eg.db', force = TRUE)
BiocManager::install('clusterProfiler', force = TRUE)
BiocManager::install("TxDb.Btaurus.UCSC.bosTau9.refGene", force = TRUE)
BiocManager::install("BSgenome.Btaurus.UCSC.bosTau9", force = TRUE)
BiocManager::install("BSgenome", force = TRUE)
BiocManager::install(version = "3.10")
BiocManager::install("Biobase", force = TRUE)
BiocManager::install(c("pathview", "gage", "gageData", "GenomicAlignments",
                       "TxDb.Hsapiens.UCSC.hg19.knownGene"), force = TRUE)
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", force = TRUE)

BiocManager::install("org.Hs.eg.db", force = TRUE)
BiocManager::install("annotate", force = TRUE)
BiocManager::install("AnnotationDbi", force = TRUE)
BiocManager::install("backports", force = TRUE)
BiocManager::install(c(
  "annotate", "AnnotationDbi", "backports", "ballgown", "BH", "Biobase", "BiocGenerics",
  "BiocManager", "BiocParallel", "BiocVersion", "Biostrings", "bit", "blob", "DBI",
  "DelayedArray", "digest", "ellipsis", "genefilter", "GenomeInfoDb", "GenomeInfoDbData",
  "GenomicAlignments", "GenomicRanges", "IRanges", "lambda.r", "limma", "pkgconfig",
  "prettyunits", "Rcpp", "RCurl", "Rhtslib", "rlang", "Rsamtools", "RSQLite",
  "rtracklayer", "S4Vectors", "SummarizedExperiment", "sva", "vctrs", "XML", "XVector",
  "zlibbioc"
), update = TRUE, ask = FALSE, force = TRUE)
BiocManager::install("methylKit")


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
getwd()
setwd("/lab_data/hansen_lab/jkincade")

date <- "06.03.2022"
folder <- "/lab_data/hansen_lab/jkincade/Methyl-Seq/Methyl-Seq_Data/DMRs/methylKit/"

bam.list <- list(
  "2132_control_R1_val_1_bismark_bt2_pe.bam",
  "2147_control_R1_val_1_bismark_bt2_pe.bam",
  "2149_control_R1_val_1_bismark_bt2_pe.bam",
  "2156_control_R1_val_1_bismark_bt2_pe.bam",
  "2163_control_R1_val_1_bismark_bt2_pe.bam",
  "2103_TI_R1_val_1_bismark_bt2_pe.bam",
  "2106_TI_R1_val_1_bismark_bt2_pe.bam",
  "2169_TI_R1_val_1_bismark_bt2_pe.bam",
  "2183_TI_R1_val_1_bismark_bt2_pe.bam",
  "2186_TI_R1_val_1_bismark_bt2_pe.bam",
  "6_146_PI_R1_val_1_bismark_bt2_pe.bam",
  "6_186_PI_R1_val_1_bismark_bt2_pe.bam",
  "7_138_PI_R1_val_1_bismark_bt2_pe.bam",
  "144_PI_R1_val_1_bismark_bt2_pe.bam",
  "198_PI_R1_val_1_bismark_bt2_pe.bam"
)

ids <- list("2132.C", "2147.C", "2149.C", "2156.C", "2163.C",
            "2103.TI", "2106.TI", "2169.TI", "2183.TI", "2186.TI",
            "6_146.PI", "6_186.PI", "7_138.PI", "144.PI", "198.PI")
treatment.list <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2) # controls = 0; TI = 1; PI = 2
class(treatment.list[1L])

# for supercomputer use - can parallel to multiple cores. 
registerDoParallel(15)

foreach(i=1:15) %dopar% {
  parallel_CpG_call<- processBismarkAln(location=bam.list[i],
                                        sample.id=ids[i],
                                        assembly = "bosTau9", save.folder = "/lab_data/hansen_lab/hanahm1/Methyl-Seq/Methyl-Seq_Data/BAM_Alignments/Called_Parallel2",
                                        save.context = "CpG", read.context = "CpG", treatment = treatment.list[i])
}

foreach(i=1:7) %dopar% {
  parallel_CHH_call<- processBismarkAln(location=bam.list[i],
                                        sample.id=ids[i],
                                        assembly = "bosTau9", save.folder = "/lab_data/hansen_lab/hanahm1/Methyl-Seq/Methyl-Seq_Data/BAM_Alignments/Called_Parallel2",
                                        save.context = "CHH", read.context = "CHH", treatment = treatment.list[i])
}

foreach(i=1:7) %dopar% {
  parallel_CHG_call<- processBismarkAln(location=bam.list[i],
                                        sample.id=ids[i],
                                        assembly = "bosTau9", save.folder = "/lab_data/hansen_lab/hanahm1/Methyl-Seq/Methyl-Seq_Data/BAM_Alignments/Called_Parallel2",
                                        save.context = "CHG", read.context = "CHG", treatment = treatment.list[i])
}

setwd("/lab_data/hansen_lab/hanahm1/Methyl-Seq/Methyl-Seq_Data/BAM_Alignments/Called_Parallel2")
CpG.Call.List <- list.files(pattern = ".*\\CpG.txt$")
CpG.Call.List <- as.list(CpG.Call.List)
CHH.Call.List <- list.files(pattern = ".*\\CHH.txt$")
CHH.Call.List <- as.list(CHH.Call.List)
CHG.Call.List <- list.files(pattern = ".*\\CHG.txt$")
CHG.Call.List <- as.list(CHG.Call.List)

CpG.Read <- methRead(CpG.Call.List, sample.id = list("272.C", "268.C", "278.C", "250.C", "842.PI", "307.PI", "308.PI"),
                     assembly = "bosTau9", treatment = c(0,0,0,0,1,1,1), context = "CpG")

CHH.Read <- methRead(CHH.Call.List, sample.id = list("272.C", "268.C", "278.C", "250.C", "842.PI", "307.PI", "308.PI"),
                     assembly = "bosTau9", treatment = c(0,0,0,0,1,1,1), context = "CHH")

CHG.Read <- methRead(CHG.Call.List, sample.id = list("272.C", "268.C", "278.C", "250.C", "842.PI", "307.PI", "308.PI"),
                     assembly = "bosTau9", treatment = c(0,0,0,0,1,1,1), context = "CHG")

getMethylationStats(CpG.Read[[2]], plot = FALSE, both.strands = FALSE)

getMethylationStats(CpG.Read[[1]], plot = TRUE, both.strands = FALSE)
getMethylationStats(CpG.Read[[2]], plot = TRUE, both.strands = FALSE)
getMethylationStats(CpG.Read[[3]], plot = TRUE, both.strands = FALSE)
getMethylationStats(CpG.Read[[4]], plot = TRUE, both.strands = FALSE)
getMethylationStats(CpG.Read[[5]], plot = TRUE, both.strands = FALSE)
getMethylationStats(CpG.Read[[6]], plot = TRUE, both.strands = FALSE)
getMethylationStats(CpG.Read[[7]], plot = TRUE, both.strands = FALSE)

getCoverageStats(CpG.Read[[2]], plot = TRUE, both.strands = FALSE)

filtered.CpG.Read <- filterByCoverage(CpG.Read, lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
filtered.CHH.Read <- filterByCoverage(CHH.Read, lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
filtered.CHG.Read <- filterByCoverage(CHG.Read, lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)

filtered.united.CpG <- unite(filtered.CpG.Read, destrand = FALSE) 
filtered.united.CHH <- unite(filtered.CHH.Read, destrand = FALSE)
filtered.united.CHG <- unite(filtered.CHG.Read, destrand = FALSE)

head(filtered.united.CpG)

getCorrelation(filtered.united.CpG, plot = TRUE)
getCorrelation(filtered.united.CHH, plot = TRUE)
getCorrelation(filtered.united.CHG, plot = TRUE)

## Filtering by standard deviations ##
percent.meth.CpG <- percMethylation(filtered.united.CpG)
percent.meth.CHH <- percMethylation(filtered.united.CHH)
percent.meth.CHG <- percMethylation(filtered.united.CHG)

CpG.matrix <- matrixStats::rowSds(percent.meth.CpG)
CHH.matrix <- matrixStats::rowSds(percent.meth.CHH)
CHG.matrix <- matrixStats::rowSds(percent.meth.CHG)

head(filtered.united.CpG[CpG.matrix > 20, ])
hist(CpG.matrix,col = 'cornflowerblue', xlab = 'Std. Dev. per CpG')

######################################################################

## Clustering Samples and PCA plots ##
clusterSamples(filtered.united.CpG, dist = "correlation", method = "ward", plot = TRUE)
clusterSamples(filtered.united.CHH, dist = "correlation", method = "ward", plot = TRUE)
clusterSamples(filtered.united.CHG, dist = "correlation", method = "ward", plot = TRUE)

PCASamples(filtered.united.CpG, screeplot = TRUE)
PCASamples(filtered.united.CpG)
PCASamples(filtered.united.CHH)
PCASamples(filtered.united.CHG)

## Calculating DMR ##
# Fisher's Exact Test #
getSampleID(filtered.united.CpG)
pooled.CpG <- pool(filtered.united.CpG, sample.ids = c("C", "PI"), treatment = c(0, 1))
CpG.Fish <- calculateDiffMeth(pooled.CpG, mc.cores = 30)

CpG.pooled.fish.all <- getMethylDiff(CpG.Fish, difference = 25, qvalue = 0.01, type = "all")
head(CpG.pooled.fish.all)

diffMethPerChr(CpG.Fish, plot = TRUE, qvalue.cutoff = 0.01, meth.cutoff = 25)

# Logistic Regression #
CpG.Diff <- calculateDiffMeth(filtered.united.CpG, mc.cores = 30)

CHH.Diff <- calculateDiffMeth(filtered.united.CHH, mc.cores = 30)

CHG.Diff <- calculateDiffMeth(filtered.united.CHG, mc.cores = 30)

diffMethPerChr(CpG.Diff, plot = TRUE, qvalue.cutoff = 0.01, meth.cutoff = 25)

CpG.hyper <- getMethylDiff(CpG.Diff, difference = 25, qvalue = 0.01, type = "hyper")
CpG.hypo <- getMethylDiff(CpG.Diff, difference = 25, qvalue = 0.01, type = "hypo")
CpG.all <-  getMethylDiff(CpG.Diff, difference = 25, qvalue = 0.01)

CHH.hyper <- getMethylDiff(CHH.Diff, difference = 25, qvalue = 0.01, type = "hyper")
CHH.hypo <- getMethylDiff(CHH.Diff, difference = 25, qvalue = 0.01, type = "hypo")
CHH.all <-  getMethylDiff(CHH.Diff, difference = 25, qvalue = 0.01)

CHG.hyper <- getMethylDiff(CHG.Diff, difference = 25, qvalue = 0.01, type = "hyper")
CHG.hypo <- getMethylDiff(CHG.Diff, difference = 25, qvalue = 0.01, type = "hypo")
CHG.all <-  getMethylDiff(CHG.Diff, difference = 25, qvalue = 0.01)
head(CpG.all)

bov <- readTranscriptFeatures("/lab_data/hansen_lab/hanahm1/bosTau9.refseq.bed.txt")

annotateWithGeneParts(as(CpG.hyper, "GRanges"),bov)
annotateWithGeneParts(as(CpG.hypo, "GRanges"),bov)
annotateWithGeneParts(as(CpG.all, "GRanges"),bov)
head(CpG.all)

annotateWithGeneParts(as(CHH.hyper, "GRanges"),bov)
annotateWithGeneParts(as(CHH.hypo, "GRanges"),bov)
annotateWithGeneParts(as(CHH.all, "GRanges"),bov)

annotateWithGeneParts(as(CHG.hyper, "GRanges"),bov)
annotateWithGeneParts(as(CHG.hypo, "GRanges"),bov)
annotateWithGeneParts(as(CHG.all, "GRanges"),bov)

annotated.CpG.hyper <- annotateWithGeneParts(as(CpG.hyper, "GRanges"), bov)
annotated.CpG.hypo <- annotateWithGeneParts(as(CpG.hypo, "GRanges"), bov)
annotated.CpG.all <- annotateWithGeneParts(as(CpG.all, "GRanges"), bov)

annotated.CHH.hyper <- annotateWithGeneParts(as(CHH.hyper, "GRanges"), bov)
annotated.CHH.hypo <- annotateWithGeneParts(as(CHH.hypo, "GRanges"), bov)
annotated.CHH.all <- annotateWithGeneParts(as(CHH.all, "GRanges"), bov)

annotated.CHG.hyper <- annotateWithGeneParts(as(CHG.hyper, "GRanges"), bov)
annotated.CHG.hypo <- annotateWithGeneParts(as(CHG.hypo, "GRanges"), bov)
annotated.CHG.all <- annotateWithGeneParts(as(CHG.all, "GRanges"), bov)

head(getAssociationWithTSS(annotated.CpG.all))
association.CpG.hyper <- getAssociationWithTSS(annotated.CpG.hyper)
class(association.CpG.all)
View(association.CpG.all)
genomation::getTargetAnnotationStats(annotated.CpG.all, percentage = TRUE, precedence = TRUE)
genomation::plotTargetAnnotation(annotated.CpG.all, precedence = TRUE, main = "CpG.All.DMR")
genomation::getFeatsWithTargetsStats(annotated.CpG.all, percentage = TRUE)

association.CpG.hyper <- getAssociationWithTSS(annotated.CpG.hyper)
association.CpG.hypo <- getAssociationWithTSS(annotated.CpG.hypo)
association.CpG.all <- getAssociationWithTSS(annotated.CpG.all)

association.CHH.hyper <- getAssociationWithTSS(annotated.CHH.hyper)
association.CHH.hypo <- getAssociationWithTSS(annotated.CHH.hypo)
association.CHH.all <- getAssociationWithTSS(annotated.CHH.all)

association.CHG.hyper <- getAssociationWithTSS(annotated.CHG.hyper)
association.CHG.hypo <- getAssociationWithTSS(annotated.CHG.hypo)
association.CHG.all <- getAssociationWithTSS(annotated.CHG.all)
dim(association.CpG.hyper)

df.CpG.hyper <- getData(CpG.hyper)
head(df.CpG.hyper)
dim(df.CpG.hyper)
#df.CpG.hyper
# for row 12, chrUN is an unknown sequence and messes with downstream analysis. Will delete
#df.CpG.hyper.test <- ifelse(grepl("chrUn",df.CpG.hyper$chr))
rowlist <- list()
L <- 1
for(i in 1:nrow(df.CpG.hyper)){
  chr <- df.CpG.hyper$chr [i]
  if(grepl("chrUn",chr)){
    rowlist [[L]] <- i
    L <- L + 1
  }
}
rowlist
rw.length <- length(rowlist)
rw.length.minus <- rowlist[[1]] - 1
rowlist[[rw.length]]
df.CpG.hyper <- df.CpG.hyper[-c(rowlist[[1]]:rowlist[[rw.length]]),]
dim(df.CpG.hyper)
df.CpG.hyper$REFSEQ <- association.CpG.hyper$feature.name
View(df.CpG.hyper)
tail(df.CpG.hyper)
df.CpG.hyper <- separate(df.CpG.hyper, "REFSEQ", c("REFSEQ", "version"), sep = "\\.")
CpG.hyper.gene.id <- bitr(df.CpG.hyper$REFSEQ, fromType = "REFSEQ", toType = c("SYMBOL", "GENENAME", "ENTREZID"), OrgDb = "org.Bt.eg.db", drop = TRUE)
symbol.CpG.hyper <- merge(CpG.hyper.gene.id, df.CpG.hyper, by = "REFSEQ")
#symbol.CpG.hyper <- merge(symbol.CpG.hyper, df.CpG.hyper.heatmap [,c(3:7)], by = "end")
head(symbol.CpG.hyper)
View(symbol.CpG.hyper)

## CpG Hypo ##
df.CpG.hypo <- getData(CpG.hypo)
head(df.CpG.hypo)
dim(df.CpG.hypo)
#df.CpG.hypo
# for row 12, chrUN is an unknown sequence and messes with downstream analysis. Will delete
#df.CpG.hypo.test <- ifelse(grepl("chrUn",df.CpG.hypo$chr))
rowlist <- list()
L <- 1
for(i in 1:nrow(df.CpG.hypo)){
  chr <- df.CpG.hypo$chr [i]
  if(grepl("chrUn",chr)){
    rowlist [[L]] <- i
    L <- L + 1
  }
}
rowlist
rw.length <- length(rowlist)
rw.length.minus <- rowlist[[1]] - 1
rowlist[[rw.length]]
df.CpG.hypo <- df.CpG.hypo[-c(rowlist[[1]]:rowlist[[rw.length]]),]
dim(df.CpG.hypo)
dim(association.CpG.hypo)
df.CpG.hypo$REFSEQ <- association.CpG.hypo$feature.name
df.CpG.hypo <- separate(df.CpG.hypo, "REFSEQ", c("REFSEQ", "version"), sep = "\\.")
CpG.hypo.gene.id <- bitr(df.CpG.hypo$REFSEQ, fromType = "REFSEQ", toType = c("SYMBOL", "GENENAME", "ENTREZID"), OrgDb = "org.Bt.eg.db", drop = TRUE)
symbol.CpG.hypo <- merge(CpG.hypo.gene.id, df.CpG.hypo, by = "REFSEQ")
#symbol.CpG.hypo <- merge(symbol.CpG.hypo, df.CpG.hypo.heatmap [,c(3:7)], by = "end")
head(symbol.CpG.hypo)


## CpG All ##
df.CpG.all <- getData(CpG.all)
head(df.CpG.all)
dim(df.CpG.all)
#df.CpG.all
# for row 12, chrUN is an unknown sequence and messes with downstream analysis. Will delete
#df.CpG.all.test <- ifelse(grepl("chrUn",df.CpG.all$chr))
rowlist <- list()
L <- 1
for(i in 1:nrow(df.CpG.all)){
  chr <- df.CpG.all$chr [i]
  if(grepl("chrUn",chr)){
    rowlist [[L]] <- i
    L <- L + 1
  }
}
rowlist
rw.length <- length(rowlist)
rw.length.minus <- rowlist[[1]] - 1
rowlist[[rw.length]]
df.CpG.all <- df.CpG.all[-c(rowlist[[1]]:rowlist[[rw.length]]),]
dim(df.CpG.all)
dim(association.CpG.all)
df.CpG.all$REFSEQ <- association.CpG.all$feature.name
df.CpG.all <- separate(df.CpG.all, "REFSEQ", c("REFSEQ", "version"), sep = "\\.")
CpG.all.gene.id <- bitr(df.CpG.all$REFSEQ, fromType = "REFSEQ", toType = c("SYMBOL", "GENENAME", "ENTREZID"), OrgDb = "org.Bt.eg.db", drop = TRUE)
symbol.CpG.all <- merge(CpG.all.gene.id, df.CpG.all, by = "REFSEQ")
#symbol.CpG.all <- merge(symbol.CpG.all, df.CpG.all.heatmap [,c(3:7)], by = "end")
head(symbol.CpG.all)

## CHG Hyper ##
df.CHG.hyper <- getData(CHG.hyper)
head(df.CHG.hyper)
dim(df.CHG.hyper)
#df.CHG.hyper
# for row 12, chrUN is an unknown sequence and messes with downstream analysis. Will delete
#df.CHG.hyper.test <- ifelse(grepl("chrUn",df.CHG.hyper$chr))
rowlist <- list()
L <- 1
for(i in 1:nrow(df.CHG.hyper)){
  chr <- df.CHG.hyper$chr [i]
  if(grepl("chrUn",chr)){
    rowlist [[L]] <- i
    L <- L + 1
  }
}
rowlist
rw.length <- length(rowlist)
rw.length.minus <- rowlist[[1]] - 1
rowlist[[rw.length]]
df.CHG.hyper <- df.CHG.hyper[-c(rowlist[[1]]:rowlist[[rw.length]]),]
dim(df.CHG.hyper)
dim(association.CHG.hyper)
df.CHG.hyper$REFSEQ <- association.CHG.hyper$feature.name
df.CHG.hyper <- separate(df.CHG.hyper, "REFSEQ", c("REFSEQ", "version"), sep = "\\.")
CHG.hyper.gene.id <- bitr(df.CHG.hyper$REFSEQ, fromType = "REFSEQ", toType = c("SYMBOL", "GENENAME", "ENTREZID"), OrgDb = "org.Bt.eg.db", drop = TRUE)
symbol.CHG.hyper <- merge(CHG.hyper.gene.id, df.CHG.hyper, by = "REFSEQ")
#symbol.CHG.hyper <- merge(symbol.CHG.hyper, df.CHG.hyper.heatmap [,c(3:7)], by = "end")
head(symbol.CHG.hyper)

## CHG Hypo ##
df.CHG.hypo <- getData(CHG.hypo)
head(df.CHG.hypo)
dim(df.CHG.hypo)
#df.CHG.hypo
# for row 12, chrUN is an unknown sequence and messes with downstream analysis. Will delete
#df.CHG.hypo.test <- ifelse(grepl("chrUn",df.CHG.hypo$chr))
rowlist <- list()
L <- 1
for(i in 1:nrow(df.CHG.hypo)){
  chr <- df.CHG.hypo$chr [i]
  if(grepl("chrUn",chr)){
    rowlist [[L]] <- i
    L <- L + 1
  }
}
rowlist
rw.length <- length(rowlist)
rw.length.minus <- rowlist[[1]] - 1
rowlist[[rw.length]]
df.CHG.hypo <- df.CHG.hypo[-c(rowlist[[1]]:rowlist[[rw.length]]),]
dim(df.CHG.hypo)
dim(association.CHG.hypo)
df.CHG.hypo$REFSEQ <- association.CHG.hypo$feature.name
df.CHG.hypo <- separate(df.CHG.hypo, "REFSEQ", c("REFSEQ", "version"), sep = "\\.")
CHG.hypo.gene.id <- bitr(df.CHG.hypo$REFSEQ, fromType = "REFSEQ", toType = c("SYMBOL", "GENENAME", "ENTREZID"), OrgDb = "org.Bt.eg.db", drop = TRUE)
symbol.CHG.hypo <- merge(CHG.hypo.gene.id, df.CHG.hypo, by = "REFSEQ")
#symbol.CHG.hypo <- merge(symbol.CHG.hypo, df.CHG.hypo.heatmap [,c(3:7)], by = "end")
head(symbol.CHG.hypo)


## CHG All ##
df.CHG.all <- getData(CHG.all)
head(df.CHG.all)
dim(df.CHG.all)
#df.CHG.all
# for row 12, chrUN is an unknown sequence and messes with downstream analysis. Will delete
#df.CHG.all.test <- ifelse(grepl("chrUn",df.CHG.all$chr))
rowlist <- list()
L <- 1
for(i in 1:nrow(df.CHG.all)){
  chr <- df.CHG.all$chr [i]
  if(grepl("chrUn",chr)){
    rowlist [[L]] <- i
    L <- L + 1
  }
}
rowlist
rw.length <- length(rowlist)
rw.length.minus <- rowlist[[1]] - 1
rowlist[[rw.length]]
df.CHG.all <- df.CHG.all[-c(rowlist[[1]]:rowlist[[rw.length]]),]
dim(df.CHG.all)
dim(association.CHG.all)
df.CHG.all$REFSEQ <- association.CHG.all$feature.name
df.CHG.all <- separate(df.CHG.all, "REFSEQ", c("REFSEQ", "version"), sep = "\\.")
CHG.all.gene.id <- bitr(df.CHG.all$REFSEQ, fromType = "REFSEQ", toType = c("SYMBOL", "GENENAME", "ENTREZID"), OrgDb = "org.Bt.eg.db", drop = TRUE)
symbol.CHG.all <- merge(CHG.all.gene.id, df.CHG.all, by = "REFSEQ")
#symbol.CHG.all <- merge(symbol.CHG.all, df.CHG.all.heatmap [,c(3:7)], by = "end")
head(symbol.CHG.all)


## CHH Hyper ##
df.CHH.hyper <- getData(CHH.hyper)
head(df.CHH.hyper)
dim(df.CHH.hyper)
#df.CHH.hyper
# for row 12, chrUN is an unknown sequence and messes with downstream analysis. Will delete
#df.CHH.hyper.test <- ifelse(grepl("chrUn",df.CHH.hyper$chr))
rowlist <- list()
L <- 1
for(i in 1:nrow(df.CHH.hyper)){
  chr <- df.CHH.hyper$chr [i]
  if(grepl("chrUn",chr)){
    rowlist [[L]] <- i
    L <- L + 1
  }
}
rowlist
rw.length <- length(rowlist)
rw.length.minus <- rowlist[[1]] - 1
rowlist[[rw.length]]
df.CHH.hyper <- df.CHH.hyper[-c(rowlist[[1]]:rowlist[[rw.length]]),]
dim(df.CHH.hyper)
dim(association.CHH.hyper)
df.CHH.hyper$REFSEQ <- association.CHH.hyper$feature.name
df.CHH.hyper <- separate(df.CHH.hyper, "REFSEQ", c("REFSEQ", "version"), sep = "\\.")
CHH.hyper.gene.id <- bitr(df.CHH.hyper$REFSEQ, fromType = "REFSEQ", toType = c("SYMBOL", "GENENAME", "ENTREZID"), OrgDb = "org.Bt.eg.db", drop = TRUE)
symbol.CHH.hyper <- merge(CHH.hyper.gene.id, df.CHH.hyper, by = "REFSEQ")
#symbol.CHH.hyper <- merge(symbol.CHH.hyper, df.CHH.hyper.heatmap [,c(3:7)], by = "end")
head(symbol.CHH.hyper)

## CHH Hypo ##
df.CHH.hypo <- getData(CHH.hypo)
head(df.CHH.hypo)
dim(df.CHH.hypo)
#df.CHH.hypo
# for row 12, chrUN is an unknown sequence and messes with downstream analysis. Will delete
#df.CHH.hypo.test <- ifelse(grepl("chrUn",df.CHH.hypo$chr))
rowlist <- list()
L <- 1
for(i in 1:nrow(df.CHH.hypo)){
  chr <- df.CHH.hypo$chr [i]
  if(grepl("chrUn",chr)){
    rowlist [[L]] <- i
    L <- L + 1
  }
}
rowlist #If empty, skip to the dim lines and make sure the row dimensions match
rw.length <- length(rowlist)
rw.length.minus <- rowlist[[1]] - 1
rowlist[[rw.length]]
df.CHH.hypo <- df.CHH.hypo[-c(rowlist[[1]]:rowlist[[rw.length]]),]
dim(df.CHH.hypo)
dim(association.CHH.hypo)
df.CHH.hypo$REFSEQ <- association.CHH.hypo$feature.name
df.CHH.hypo <- separate(df.CHH.hypo, "REFSEQ", c("REFSEQ", "version"), sep = "\\.")
CHH.hypo.gene.id <- bitr(df.CHH.hypo$REFSEQ, fromType = "REFSEQ", toType = c("SYMBOL", "GENENAME", "ENTREZID"), OrgDb = "org.Bt.eg.db", drop = TRUE)
symbol.CHH.hypo <- merge(CHH.hypo.gene.id, df.CHH.hypo, by = "REFSEQ")
#symbol.CHH.hypo <- merge(symbol.CHH.hypo, df.CHH.hypo.heatmap [,c(3:7)], by = "end")
head(symbol.CHH.hypo)


## CHH All ##
df.CHH.all <- getData(CHH.all)
head(df.CHH.all)
dim(df.CHH.all)
#df.CHH.all
# for row 12, chrUN is an unknown sequence and messes with downstream analysis. Will delete
#df.CHH.all.test <- ifelse(grepl("chrUn",df.CHH.all$chr))
rowlist <- list()
L <- 1
for(i in 1:nrow(df.CHH.all)){
  chr <- df.CHH.all$chr [i]
  if(grepl("chrUn",chr)){
    rowlist [[L]] <- i
    L <- L + 1
  }
}
rowlist
rw.length <- length(rowlist)
rw.length.minus <- rowlist[[1]] - 1
rowlist[[rw.length]]
df.CHH.all <- df.CHH.all[-c(rowlist[[1]]:rowlist[[rw.length]]),]
dim(df.CHH.all)
dim(association.CHH.all)
df.CHH.all$REFSEQ <- association.CHH.all$feature.name
df.CHH.all <- separate(df.CHH.all, "REFSEQ", c("REFSEQ", "version"), sep = "\\.")
CHH.all.gene.id <- bitr(df.CHH.all$REFSEQ, fromType = "REFSEQ", toType = c("SYMBOL", "GENENAME", "ENTREZID"), OrgDb = "org.Bt.eg.db", drop = TRUE)
symbol.CHH.all <- merge(CHH.all.gene.id, df.CHH.all, by = "REFSEQ")
#symbol.CHH.all <- merge(symbol.CHH.all, df.CHH.all.heatmap [,c(3:7)], by = "end")
head(symbol.CHH.all)

symbol.CpG.hyper$context <- "CpG"
symbol.CpG.hyper$state <- "Hyper"
symbol.CpG.hypo$context <- "CpG"
symbol.CpG.hypo$state <- "Hypo"

symbol.CHG.hyper$context <- "CHG"
symbol.CHG.hyper$state <- "Hyper"
symbol.CHG.hypo$context <- "CHG"
symbol.CHG.hypo$state <- "Hypo"

symbol.CHH.hyper$context <- "CHH"
symbol.CHH.hyper$state <- "Hyper"
symbol.CHH.hypo$context <- "CHH"
symbol.CHH.hypo$state <- "Hypo"

merged <- rbind(symbol.CpG.hyper, symbol.CpG.hypo)
merged <- rbind(merged, symbol.CHG.hyper)
merged <- rbind(merged, symbol.CHG.hypo)
merged <- rbind(merged, symbol.CHH.hyper)
merged <- rbind(merged, symbol.CHH.hypo)

merged$Pred.Expr <- merged$meth.diff * -1

write.csv(merged, file = paste(folder, "merged-DMRs", date, ".csv", sep = ""), row.names = FALSE)


write.csv(symbol.CpG.hypo, file = paste(folder, "CpG-Hypo_", date, ".csv", sep = ""), row.names = FALSE)
write.csv(symbol.CpG.hyper, file = paste(folder, "CpG-hyper_", date, ".csv", sep = ""), row.names = FALSE)
write.csv(symbol.CpG.all, file = paste(folder, "CpG-all_", date, ".csv", sep = ""), row.names = FALSE)

write.csv(symbol.CHG.hypo, file = paste(folder, "CHG-hypo_", date, ".csv", sep = ""), row.names = FALSE)
write.csv(symbol.CHG.hyper, file = paste(folder, "CHG-hyper_", date, ".csv", sep = ""), row.names = FALSE)
write.csv(symbol.CHG.all, file = paste(folder, "CHG-all_", date, ".csv", sep = ""), row.names = FALSE)

write.csv(symbol.CHH.hypo, file = paste(folder, "CHH-hypo_", date, ".csv", sep = ""), row.names = FALSE)
write.csv(symbol.CHH.hyper, file = paste(folder, "CHH-hyper_", date, ".csv", sep = ""), row.names = FALSE)
write.csv(symbol.CHH.all, file = paste(folder, "CHH-all_", date, ".csv", sep = ""), row.names = FALSE)

merged.CpG <- cbind(symbol.CpG.hypo, symbol.CpG.hyper)

View(symbol.CpG.hypo)

#####################################################################################
## GO and KEGG ##

CpG.all.GO <- groupGO(gene = symbol.CpG.all$ENTREZID, "org.Bt.eg.db", ont = "BP", level = 3, readable = TRUE)
CpG.hyper.GO <- groupGO(gene = symbol.CpG.hyper$ENTREZID, "org.Bt.eg.db", ont = "BP", level = 3, readable = TRUE)
CpG.hypo.GO <- groupGO(gene = symbol.CpG.hypo$ENTREZID, "org.Bt.eg.db", ont = "BP", level = 3, readable = TRUE)

CHH.all.GO <- groupGO(gene = symbol.CHH.all$ENTREZID, "org.Bt.eg.db", ont = "BP", level = 3, readable = TRUE)
CHH.hyper.GO <- groupGO(gene = symbol.CHH.hyper$ENTREZID, "org.Bt.eg.db", ont = "BP", level = 3, readable = TRUE)
CHH.hypo.GO <- groupGO(gene = symbol.CHH.hypo$ENTREZID, "org.Bt.eg.db", ont = "BP", level = 3, readable = TRUE)

CHG.all.GO <- groupGO(gene = symbol.CHG.all$ENTREZID, "org.Bt.eg.db", ont = "BP", level = 3, readable = TRUE)
CHG.hyper.GO <- groupGO(gene = symbol.CHG.hyper$ENTREZID, "org.Bt.eg.db", ont = "BP", level = 3, readable = TRUE)
CHG.hypo.GO <- groupGO(gene = symbol.CHG.hypo$ENTREZID, "org.Bt.eg.db", ont = "BP", level = 3, readable = TRUE)



barplot(CpG.all.GO, drop = TRUE, showCategory = 10, xlab = "Gene Count", title = "CpG.all")
barplot(CpG.hyper.GO, drop = TRUE, showCategory = 10, xlab = "Gene Count", title = "CpG.hyper")
barplot(CpG.hypo.GO, drop = TRUE, showCategory = 10, xlab = "Gene Count", title = "CpG.hypo")

barplot(CHH.all.GO, drop = TRUE, showCategory = 10, xlab = "Gene Count", title = "CHH.all")
barplot(CHH.hyper.GO, drop = TRUE, showCategory = 10, xlab = "Gene Count", title = "CHH.hyper")
barplot(CHH.hypo.GO, drop = TRUE, showCategory = 10, xlab = "Gene Count", title = "CHH.hypo")

barplot(CHG.all.GO, drop = TRUE, showCategory = 10, xlab = "Gene Count", title = "CHG.all")
barplot(CHG.hyper.GO, drop = TRUE, showCategory = 10, xlab = "Gene Count", title = "CHG.hyper")
barplot(CHG.hypo.GO, drop = TRUE, showCategory = 10, xlab = "Gene Count", title = "CHG.hypo")





#####################################################################################
# for ViewBS MethHeatmap function- feed output into bash command #
CHG.hyper.heatmap <- df.CHG.hyper.heatmap[,c(1:3)]
CHG.hyper.heatmap <-CHG.hyper.heatmap [-12,]
#colnames(CHG.hyper.heatmap) <- NULL
#CHG.hyper.heatmap3 <-CHG.hyper.heatmap [,-4]

head(symbol.CHG.hyper)
df.CHG.hyper.heatmap <- symbol.CHG.hyper[,c(4,5,6,2)]
head(df.CHG.hyper.heatmap)


CHG.hyper.heatmap$feature.id [c(1:11)] <- association.CHG.hyper$feature.name [c(1:11)]
View(CHG.hyper.heatmap)
CHG.hyper.heatmap <- separate(CHG.hyper.heatmap, "feature.id", c("feature", "version"), sep = "\\.") #this is needed for org.db to read refseq numbers
colnames(CHG.hyper.heatmap.test) [4] <- "REFSEQ"
head(CHG.hyper.heatmap.test)
#CHG.hyper.heatmap$feature.id [10] <- "XM_024986389.1.a"
#CHG.hyper.heatmap$feature.id [11] <- "XM_024986389.1.b"

write.table(df.CpG.hyper, file = "/lab_data/hansen_lab/hanahm1/Methyl-Seq/Methyl-Seq_Data/DMRs/CpG_hyper.txt", append = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

write.table(df.CHG.hyper.heatmap, file = "/lab_data/hansen_lab/hanahm1/Methyl-Seq/Methyl-Seq_Data/DMRs/CHG_hyper.txt", append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(CHG.hyper.heatmap3, file = "/lab_data/hansen_lab/hanahm1/Methyl-Seq/Methyl-Seq_Data/DMRs/CHG_hyper3.txt", append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

gene.id <- bitr(CHG.hyper.heatmap.test$feature, fromType = "REFSEQ", toType = c("SYMBOL", "GENENAME", "ENTREZID"), OrgDb = "org.Bt.eg.db", drop = TRUE)

gene.id

symbol.CHG.hyper.heatmap <- merge(gene.id, CHG.hyper.heatmap.test, by = "REFSEQ")
symbol.CHG.hyper.heatmap <- merge(symbol.CHG.hyper.heatmap, df.CHG.hyper.heatmap [,c(3:7)], by = "end")
head(symbol.CHG.hyper.heatmap)

## For when BITR doesn't work ##
row.names(CHG.hyper.heatmap) <- NULL
x <- org.Bt.egREFSEQ2EG
mapped_genes <- mappedkeys(x)
class(mapped_genes)
xx.df <- as.data.frame(x[mapped_genes])
head(xx.df)
xx<- as.list(x[mapped_genes])
head(xx)

for(i in 1:nrow(CHG.hyper.heatmap)){
  id <- CHG.hyper.heatmap$feature [i]
  for(x in 1:nrow(xx.df)){
    acc <- xx.df$accession [x]
    en.id <- xx.df$gene_id [x]
    if (id == acc){
      CHG.hyper.heatmap$ensemble [i] <- en.id
    } else{print("nope")}
  }
}

CHG.hyper.heatmap.test <- CHG.hyper.heatmap[,-6]

strsplit(CHG.hyper.heatmap$feature.id[1], ".")




