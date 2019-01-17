### Perform various RNAseq tasks: QC and DEG analysis ###
# Arg 1: folder on which to run
# Arg 2: path to GTF of genome
# To install bioconductor packages for R3.4, use source("https://bioconductor.org/biocLite.R")
#biocLite("apeglm")
#biocLite("edgeR")
#biocLite("DESeq2")
#biocLite("ReportingTools")
#biocLite("regionReport")

library(QoRTs)
library(DESeq2)
library(edgeR)
library(ReportingTools)
library(regionReport)
library(ggplot2) ## For theme_bw() and plotting in reports

args<-commandArgs(TRUE)

setwd(args[1])
QCdirs <- dir()

#### RUN QoRTs FOR BY SAMPLE QC ####
# Hard-coded sample decoder.
decoder.data <- data.frame(unique.ID = 1:length(QCdirs),
								group.ID = c("CONTROL","CONTROL2","CONTROL3"),
								sample.ID = c("20_animals_003","10_animals_001","20_animals_002"),
								qc.data.dir = QCdirs)

# Read in results from QC folders
	
res <- read.qc.results.data("./", decoder = decoder.data, calc.DESeq2 = TRUE, calc.edgeR = TRUE)

# Make multiplot, defaults
makeMultiPlot.all(res, outfile.dir = "./QoRTS_summary_", plot.device.name = "pdf")

# Make multiplot, colored by sample
makeMultiPlot.colorBySample(res, plot.device.name="pdf")


#### RUN DESeq2 ####

# Create DESeq2 Object
sampleFiles <- paste0(decoder.data$qc.data.dir,"/QC.geneCounts.formatted.for.DESeq.txt.gz")
sampleCondition <- decoder.data$group.ID
sampleName <- decoder.data$sample.ID

sampleTable <- data.frame(sampleName = sampleName,
							fileName = sampleFiles,
							condition = sampleCondition)

# Make counts table
dds_table <-  DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design = ~ condition)
dds <- DESeq(dds_table)
resultsNames(dds)

# Default DESeq2 tests
DESeqResults <- results(dds)
resOrdered <- DESeqResults[order(DESeqResults$pvalue),]
summary(resOrdered)
sum(resOrdered$padj < 0.1, na.rm=TRUE)

# Needs to be optimized for particular experimental design
DESeqResults_ControlVs2 <- results(dds, contrast=c("condition","CONTROL","CONTROL2"))
DESeqResults_ControlVs3 <- results(dds, contrast=c("condition","CONTROL","CONTROL3"))

resLFC <- lfcShrink(dds, coef="condition_CONTROL3_vs_CONTROL", type="apeglm")

pdf(file="plotMA.pdf")
plotMA(resLFC, ylim=c(-2,2), xlim=c(0,500))
dev.off()


#### Create Report for DESeq2 Results ####
## makeDESeqDF? makeDESeqDF(object, countTable, pvalueCutoff, conditions, annotation.db, expName, reportDir, ...)

Folsomia_gtf <- GenomicFeatures::makeTxDbFromGFF(args[2])

# desReport <- HTMLReport(shortName ='Folsomia_RNAseq_DESeq_test',
#						title = 'Folsomia candida RNA-seq analysis of differential expression using DESeq',
#						reportDirectory = "./reports")

## P value cutoff should really be 0.05 or lower, this is for demonstration purposes
# publish(DESeqResults, desReport, pvalueCutoff=1,
#		annotation.db=Folsomia_gtf, factor = colData(DESeqResults)$condition,
#		reportDir="./reports")

# finish(desReport)


#### Create regionReport of DESeq2 results ####

report <- DESeq2Report(dds, project = 'DESeq2 HTML report',
    intgroup = c('condition'), outdir = 'DESeq2-example',
    output = 'index', theme = theme_bw())


