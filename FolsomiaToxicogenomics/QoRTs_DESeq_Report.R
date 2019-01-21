### Perform various RNAseq tasks: QC and DEG analysis ###
# Arg 1: folder on which to run
# Arg 2: path to GTF of genome
# To install bioconductor packages for R3.4, use source("https://bioconductor.org/biocLite.R")
#biocLite("apeglm")
#biocLite("edgeR")
#biocLite("DESeq2")
#biocLite("ReportingTools")
#biocLite("regionReport")
#biocLite("clusterProfiler")
#biocLite("GOexpress")
library(GOexpress)
library(biomaRt)
library(QoRTs)
library(DESeq2)
library(edgeR)
library(ReportingTools)
library(regionReport)
library(ggplot2) ## For theme_bw() and plotting in reports

args<-commandArgs(TRUE)

setwd(args[1])
QCdirs <- dir(pattern="R_")

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

# Plot biotypes
byLane.plotter <- build.plotter.colorBySample(res)
pdf(file = "./biotype.pdf")
makePlot.biotype.rates(byLane.plotter)
dev.off()

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

# Make for GOexpress
exprs <- counts(dds, normalized=T) # get normalized counts for each feature per sample
dimnames(exprs) = list(rownames(exprs), colnames(exprs))  # set the sample names correctly

decoder.data.goexpress <- decoder.data

row.names(decoder.data.goexpress) <- colnames(exprs)
phenoData <- new("AnnotatedDataFrame", data=decoder.data.goexpress)
minimalSet <- ExpressionSet(assayData=exprs, phenoData=phenoData)


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
#		title = 'Folsomia candida RNA-seq analysis of differential expression using DESeq',
#		reportDirectory = "./reports")

## P value cutoff should really be 0.05 or lower, this is for demonstration purposes
# publish(DESeqResults, desReport, pvalueCutoff=1,
#		annotation.db=Folsomia_gtf, factor = colData(DESeqResults)$condition,
#		reportDir="./reports")

# finish(desReport)

#### Create regionReport of DESeq2 results ####

report <- DESeq2Report(dds, project = 'DESeq2 HTML report',
    intgroup = c('condition'), outdir = 'DESeq2-report',
    output = 'index', theme = theme_bw())
### BUILD ORG PACKAGE OF ANNOTATIONS FOR FOLSOMIA - SHOULD ONLY NEED TO BE DONE ONCE.
library(AnnotationForge)

# makeOrgPackageFromNCBI(version = "0.1",
#                       author = "Matt Meier <matthew.meier@canada.ca>",
#                       maintainer = "Matt Meier <matthew.meier@canada.ca>",
#                       outputDir = ".",
#                       tax_id = "158441",
#                       genus = "Folsomia",
#                       species = "candida")

### Look at gene lists for enriched pathways

# install.packages("./org.Fcandida.eg.db", repos=NULL)

library(org.Fcandida.eg.db)

### Create geneList object

d <- data.frame(row.names(DESeqResults),DESeqResults[,2])
## 1st column is gene ID
## 2nd column is FC

## feature 1: numeric vector
geneList = d[,2]
## feature 2: named vector
names(geneList) = as.character(d[,1])
## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)

genesymbols=d[,1]


# ClusterProfiler
library(OrganismDbi)
library(clusterProfiler)
library(biomaRt)

listMarts(host="metazoa.ensembl.org")
ensembl=useMart("metazoa_mart", host="metazoa.ensembl.org")
listDatasets(ensembl)
ensembl = useDataset("fcandida_eg_gene", mart=ensembl)
ensembl = useDataset("fcandida_eg_gene", mart=useMart("metazoa_mart", host="metazoa.ensembl.org"))
filters = listFilters(ensembl)
filters[1:5,]
attributes = listAttributes(ensembl)
attributes[1:5,]

# Fcandida <- makeTxDbFromBiomart(biomart="metazoa_mart", dataset="fcandida_eg_gene",  host="metazoa.ensembl.org")
# saveDb(Fcandida, file="~/dbs/folsomia/GCF_002217175.1_ASM221717v1/ensembl/Fcandida.sqlite")
Fcandida <- loadDb("Fcandida.sqlite")

for (i in keytypes(ensembl)) { print(i); print(head(keys(ensembl, keytype=i))) }

goIDs <- getBM(attributes=c('ensembl_gene_id', 'go_id'), 
      filters = 'ensembl_gene_id', 
      values = genesymbols, 
      mart = ensembl)

allgenes.Ensembl = getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id', 'description'), mart=ensembl)
colnames(allgenes.Ensembl)[1] = 'gene_id'
allGO.Ensembl = getBM(attributes=c('go_id', 'name_1006', 'namespace_1003'), mart=ensembl)
GOgenes.Ensembl = getBM(attributes=c('ensembl_gene_id', 'go_id'), mart=ensembl)
colnames(GOgenes.Ensembl)[1] = 'gene_id'
GOgenes.Ensembl = GOgenes.Ensembl[GOgenes.Ensembl$go_id != '',]
GOgenes.Ensembl = GOgenes.Ensembl[GOgenes.Ensembl$gene_id != '',]
save(GOgenes.Ensembl, file='GOgenes.Ensembl.rda')
save(allGO.Ensembl, file='allGO.Ensembl.rda')
save(allgenes.Ensembl, file='allgenes.Ensembl.rda')

go_express_analysis <- GO_analyse(eSet=minimalSet, f='sample.ID', GO_genes=GOgenes.Ensembl, all_GO=allGO.Ensembl, all_genes=allgenes.Ensembl)

BP.5 <- subset_scores(result = AlvMac_results.pVal, namespace = "biological_process", total = 5, p.val=0.05) # requires 5 or more associated genes

<- subset( go_express_analysis.pVal$GO, total_count >= 5 & p.val<0.05 & namespace_1003=='biological_process')

heatmap_GO(go_id = "GO:0016507", result = go_express_analysis, eSet=minimalSet, cexRow=0.4, cexCol=1, cex.main=1, main.Lsplit=30)
expression_plot(gene_id = "Fcan01_154", result = go_express_analysis , eSet=minimalSet, x_var = "sample.ID", title.size=1.5, legend.title.size=10, legend.text.size=10, legend.key.size=15)

# groupGO(), enrichGO(), and gseGO() may be useful here

> columns(org.Fcandida.eg.db)
[1] "ACCNUM"   "ALIAS"    "CHR"      "ENTREZID" "GENENAME" "GID"      "PMID"
[8] "REFSEQ"   "SYMBOL"
> keytypes(org.Fcandida.eg.db)
[1] "ACCNUM"   "ALIAS"    "ENTREZID" "GENENAME" "GID"      "PMID"     "REFSEQ"
[8] "SYMBOL"

for (i in keytypes(org.Fcandida.eg.db)) { print(i); print(head(keys(org.Fcandida.eg.db, keytype=i))) }
for (i in keytypes(Fcandida)) { print(i); print(head(keys(Fcandida, keytype=i))) }


[1] "ACCNUM"
[1] "AAK77866.1" "AAL67688.1" "AAL78089.1" "AB190297.1" "AB190298.1"
[6] "AB435158.1"

[1] "ALIAS"
[1] "ATG8a"      "Actin-3"    "Actin-5C"   "Adseverin"  "Afadin"
[6] "Aftiphilin"

[1] "ENTREZID"
[1] "110841585" "110841586" "110841587" "110841588" "110841589" "110841590"

[1] "GENENAME"
[1] "(E3-independent) E2 ubiquitin-conjugating enzyme UBE2O-like"
[2] "1,2-dihydroxy-3-keto-5-methylthiopentene dioxygenase-like"
[3] "1,25-dihydroxyvitamin D(3) 24-hydroxylase, mitochondrial-like"
[4] "1,4-alpha-glucan-branching enzyme-like"
[5] "1,5-anhydro-D-fructose reductase-like"
[6] "1-acyl-sn-glycerol-3-phosphate acyltransferase alpha-like"

[1] "GID"
[1] "110841586" "110841587" "110841590" "110841591" "110841603" "110841607"

[1] "PMID"
[1] "16174032" "18853822" "20022415" "20353601" "20921393" "23873479"

[1] "REFSEQ"
[1] "XM_022087273.1" "XM_022087274.1" "XM_022087275.1" "XM_022087276.1"
[5] "XM_022087277.1" "XM_022087278.1"

[1] "SYMBOL"
[1] "LOC110841585" "LOC110841586" "LOC110841587" "LOC110841588" "LOC110841589"
[6] "LOC110841590"

head(bitr(d[,1], fromType="SYMBOL", toType=c("ACCNUM"), OrgDb="org.Fcandida.eg.db"))
head(bitr(d[,1], fromType="SYMBOL", toType=c("UNIPROT"), OrgDb="org.Fcandida.eg.db")) ### Won't work because no table exists for this
### 'toType' should be one of ACCNUM, ALIAS, ENTREZID, GENENAME, GID, PMID, REFSEQ, SYMBOL.


> head(keys(org.Fcandida.eg.db, keytype="ENTREZID"))
[1] "110841585" "110841586" "110841587" "110841588" "110841589" "110841590"
> head(keys(org.Fcandida.eg.db, keytype="GID"))
[1] "110841586" "110841587" "110841590" "110841591" "110841603" "110841607"
> head(keys(org.Fcandida.eg.db, keytype="REFSEQ"))
[1] "XM_022087273.1" "XM_022087274.1" "XM_022087275.1" "XM_022087276.1"
[5] "XM_022087277.1" "XM_022087278.1"
> head(keys(org.Fcandida.eg.db, keytype="ACCNUM"))
[1] "AAK77866.1" "AAL67688.1" "AAL78089.1" "AB190297.1" "AB190298.1"
[6] "AB435158.1"
> head(keys(org.Fcandida.eg.db, keytype="SYMBOL"))
[1] "LOC110841585" "LOC110841586" "LOC110841587" "LOC110841588" "LOC110841589"
[6] "LOC110841590"

## I think the column d[,1]  needs to be converted to ENTREZID - look into this.



ggo <- groupGO(gene     = d[,1],
               OrgDb    = Fcandida,
			   keytype  = "GENEID",
               ont      = "BP",
               level    = 1,
               readable = TRUE)

			   
ggo <- groupGO(gene     = 110841585,
               OrgDb    = org.Fcandida.eg.db,
               ont      = "BP",
               level    = 1,
               readable = TRUE)

head(ggo)

ego <- gseGO(geneList     = geneList,
              OrgDb        = org.Fcandida.eg.db,
              ont          = "CC",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
			  
barplot(ego, drop=TRUE, showCategory=12)
barplot(ggo, drop=TRUE, showCategory=12)
