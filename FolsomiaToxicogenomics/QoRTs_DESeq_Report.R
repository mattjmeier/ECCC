### Perform various RNAseq tasks: QC and DEG analysis ###
# Arg 1: folder on which to run
# Arg 2: path to GTF of genome

# To install bioconductor packages for R3.4...
#source("https://bioconductor.org/biocLite.R")
#biocLite("apeglm")
#biocLite("edgeR")
#biocLite("DESeq2")
#biocLite("ReportingTools")
#biocLite("regionReport")
#biocLite("clusterProfiler")
#biocLite("GOexpress")
#biocLite("QoRTs")
#biocLite("RMariaDB")
#biocLite("vsn")

## Load required libraries...
# library(OrganismDbi) ## Only load if making a new organism annotation package
# library(AnnotationForge) ## Only load if making a new organism annotation package
# library(GOexpress)
# library(ReportingTools)

packages <- c("RMariaDB",
                "clusterProfiler",
                "biomaRt",
                "QoRTs",
                "DESeq2",
                "edgeR",
                "regionReport",
                "org.Fcandida.eg.db",
                "ggplot2",
                "magrittr",
                "vsn"
                )

invisible(suppressPackageStartupMessages(lapply(packages, function(x)require(x, character.only = T, quietly = T))))

args<-commandArgs(TRUE)

setwd(args[1]) # QC folder
QCdirs <- dir()

###########################################################################################
####################################
#### RUN QoRTs FOR BY-SAMPLE QC ####
####################################
###########################################################################################

# Hard-coded sample decoder. Should ideally be provided as a text file to be read in.
decoder.data <- read.table("./metadata.ordered.txt", header=T, sep="\t") 

# Read in results from QC folders

res <- read.qc.results.data("./", decoder = decoder.data, calc.DESeq2 = TRUE, calc.edgeR = TRUE)

# Make multiplot, defaults
#makeMultiPlot.all(res, outfile.dir = "./QoRTS_summary_", plot.device.name = "pdf")

# Make multiplot, colored by sample
#makeMultiPlot.colorBySample(res, plot.device.name="pdf")

# Plot biotypes
byLane.plotter <- build.plotter.colorBySample(res)
#pdf(file = "./biotype.pdf")
#makePlot.biotype.rates(byLane.plotter)
#dev.off()

###########################################################################################
###########################################################################################
#### RUN DESeq2 ###########################################################################
###########################################################################################
###########################################################################################

# Create DESeq2 Object
sampleFiles <- paste0(decoder.data$qc.data.dir,"/QC.geneCounts.formatted.for.DESeq.txt.gz")
sampleCondition <- as.factor(decoder.data$group.ID)
sampleName <- decoder.data$sample.ID

sampleTable <- data.frame(sampleName = sampleName,
                          fileName = sampleFiles,
                          condition = sampleCondition)

# Create DESeq object
dds_from_counts <-  DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design = ~ condition)
### Runs DGE analysis ###
dds <- DESeq(dds_from_counts)
# Default DESeq2 tests
DESeqResults <- results(dds)
resOrdered <- DESeqResults[order(DESeqResults$pvalue),]
summary(resOrdered)
sum(resOrdered$padj < 0.1, na.rm=TRUE)

# Needs to be customized for particular experimental design
#DESeqResults_ControlVs2 <- results(dds, contrast=c("condition","CONTROL","CONTROL2"))
#DESeqResults_ControlVs3 <- results(dds, contrast=c("condition","CONTROL","CONTROL3"))
# Logfoldchange Shrink
#resLFC <- lfcShrink(dds, coef="condition_CONTROL3_vs_CONTROL", type="apeglm")

# Some info for sanity check
resultsNames(dds)
colData(dds) %>% head
assay(dds) %>% head
rowRanges(dds) %>% head
counts(dds) %>% str

# Make read counts table
read.counts <- counts(dds, normalized=F)
# Normalize read counts by size factors
read.counts.sf_normalized <- counts(dds, normalized=T)
lognorm.read.counts <- log2(read.counts.sf_normalized + 1)
# Rlog normalized
dds.rlog <- rlog(dds, blind=T) ### May need to set blind=F if large global changes are observed
rlog.norm.counts <- assay(dds.rlog)
# Calculate distances, Pearson correlation
distance.m_rlog <- as.dist(1 - cor(rlog.norm.counts, method="pearson"))
cor(rlog.norm.counts, method="pearson")
# Principal components analysis
pc <- prcomp(t(rlog.norm.counts))

###########################################################################################
###########################################################################################
### RUN edgeR #############################################################################
###########################################################################################
###########################################################################################
### IMPORTANT NOTE: ### This will ONLY work with an appropriately designed experiment. ####
####################### You must have replicates in each group to run these tests. ########
###########################################################################################

# read.counts matrix as input...
# sampleCondition required for grouping...

group <- factor(sampleCondition) # Just to make sure it is a factor
design <- model.matrix(~group)
y <- DGEList(counts=read.counts, group=group)
        head(y$counts)
        head(y$samples)
        summary(log2(rowSums(cpm(y))))
        # Filter based on coverage
        keep <- rowSums( cpm(y) >= 1) >= 2 # 1 CPM in at least 2 samples... Adjust accordingly
        y.keep <- y[keep,]
        summary(log2(rowSums(cpm(y.keep))))
        # Recompute library sizes
        y.keep$samples$lib.size <- colSums(y.keep$counts)
y.keep <- calcNormFactors(y.keep)
y.keep$samples
#y.keep <- estimateGLMCommonDisp(y.keep,design)
#y.keep <- estimateGLMTrendedDisp(y.keep,design)
#y.keep <- estimateGLMTagwiseDisp(y.keep,design)
#fit <- glmFit(y.keep,design)
#lrt <- glmLRT(fit,coef=2) ### coef depends on experiment
#topTags(lrt)
#edgeR_results <- topTags(lrt, n=Inf, sort.by="Pvalue", adjust.method="BH")

###########################################################################################
######## PRODUCE VARIOUS PLOTS IN PDF #####################################################
###########################################################################################

# Plot read counts and log normalized read counts for QC check
pdf(file="Read_counts.log2normalized.vs.untransformed.pdf")
par(mfrow=c(2,1))
boxplot(read.counts.sf_normalized,
            notch=T,
            main="Untransformed read counts",
            ylab="Read counts")
boxplot(lognorm.read.counts,
            notch=T,
            main="log2-Transformed read counts",
            ylab="log2 of Read counts")
par(mfrow=c(1,1))

# Plot histogram of counts per million, edgeR data
hist(log2(rowSums(cpm(y))), main="edgeR counts per million, before filtering")
hist(log2(rowSums(cpm(y.keep))), main="edgeR counts per million, after filtering")

# Plot first and second samples against each other
plot(lognorm.read.counts[,1:2],
            cex=0.5,
            main="Size factor and log normalized read counts")
plot(rlog.norm.counts[,1:2],
            cex=0.5,
            main="Size factor and rlog normalized read counts")

# Plot mean by standard deviation
msd_plot <- meanSdPlot(lognorm.read.counts, ranks=FALSE, plot=FALSE)
msd_plot$gg +
            ggtitle("Sequencing depth normalized log2(read counts)") +
            ylab("Standard deviation")
msd_plot2 <- meanSdPlot(rlog.norm.counts, ranks=FALSE, plot=FALSE)
msd_plot2$gg +
            ggtitle("Sequencing depth normalized rlog(read counts)") +
            ylab("Standard deviation")
plot(hclust(distance.m_rlog),
            labels=colnames(rlog.norm.counts),
            main="rlog transformed read counts\ndistance: Pearson correlation")

# Principal components analysis plot
P <- plotPCA(dds.rlog)
P <- P + theme_bw() + ggtitle("Rlog transformed counts")
print(P)

# P value distribution
hist(DESeqResults$pvalue, col="grey",
     border="white", xlab="", ylab="",
     main="Frequencies of all unadjusted p-values")

# MA Plots
DESeq2::plotMA(DESeqResults, alpha=0.05, main="MA plot,  results", ylim=c(-4,4))
#DESeq2::plotMA(resLFC, alpha=0.05, main="MA plot, LFC shrunk", ylim=c(-4,4))

dev.off()
###########################################################################################

# Make objects for later use in GOexpress
dimnames(read.counts.sf_normalized) = list(rownames(read.counts.sf_normalized), colnames(read.counts.sf_normalized))  # set the sample names correctly
# New metadata table for GOexpress
decoder.data.goexpress <- decoder.data
# New objects for input to GOexpress
row.names(decoder.data.goexpress) <- colnames(read.counts.sf_normalized)
phenoData <- new("AnnotatedDataFrame", data=decoder.data.goexpress)
minimalSet <- ExpressionSet(assayData=read.counts.sf_normalized, phenoData=phenoData)
#### see https://support.bioconductor.org/p/73347/


###########################################################################################
#### Create Region Report for DESeq2 and edgeR Results ####################################
###########################################################################################

# Maybe useful: makeDESeqDF(object, countTable, pvalueCutoff, conditions, annotation.db, expName, reportDir, ...)

# Folsomia_gtf <- GenomicFeatures::makeTxDbFromGFF(args[2])

# desReport <- HTMLReport(shortName ='Folsomia_RNAseq_DESeq_test',
#               title = 'Folsomia candida RNA-seq analysis of differential expression using DESeq',
#               reportDirectory = "./reports")

## P value cutoff should really be 0.05 or lower, this is for demonstration purposes
# publish(DESeqResults, desReport, pvalueCutoff=1,
#               annotation.db=Folsomia_gtf, factor = colData(DESeqResults)$condition,
#               reportDir="./reports")

# finish(desReport)

plotMA <- DESeq2::plotMA

report <- DESeq2Report(dds,
                       project = 'DESeq2 HTML report',
                       intgroup = c('condition'),
                       outdir = 'DESeq2-report',
                       output = 'index',
                       theme = theme_bw())



# edgeReport()

###########################################################################################
### BUILD ORG PACKAGE OF ANNOTATIONS FOR FOLSOMIA - SHOULD ONLY NEED TO BE DONE ONCE.######
### THERE ARE MANY WAYS TO SKIN THIS CAT. CUSTOM DATA FRAMES GIVE THE MOST FLEXIBILITY.####
### EXAMPLES FOLLOW. ######################################################################
###########################################################################################

# makeOrgPackageFromNCBI(version = "0.1",
#                       author = "Matt Meier <matthew.meier@canada.ca>",
#                       maintainer = "Matt Meier <matthew.meier@canada.ca>",
#                       outputDir = ".",
#                       tax_id = "158441",
#                       genus = "Folsomia",
#                       species = "candida")

#### LOAD FOLSOMIA BIOMART ACCESS #### Useful to do in script routinely, for possible use later

listMarts(host="metazoa.ensembl.org")
ensembl=useMart("metazoa_mart", host="metazoa.ensembl.org")
listDatasets(ensembl)
ensembl = useDataset("fcandida_eg_gene", mart=ensembl)
ensembl = useDataset("fcandida_eg_gene", mart=useMart("metazoa_mart", host="metazoa.ensembl.org"))
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
for (i in keytypes(ensembl)) { print(i); print(head(keys(ensembl, keytype=i))) }

###########################################################################################
# This works, but does not include GO IDs
# Fcandida <- makeTxDbFromBiomart(biomart="metazoa_mart",
#                                 dataset="fcandida_eg_gene",
#                                 host="metazoa.ensembl.org")
# saveDb(Fcandida, file="~/dbs/folsomia/GCF_002217175.1_ASM221717v1/ensembl/Fcandida.sqlite")
# install.packages("./TxDb.Fcandida.BioMart.metazoamart", repos = NULL, type="source")
# library(TxDb.Fcandida.BioMart.metazoamart)
###########################################################################################

###########################################################################################
# This one didn't work... Probably designed for standard Ensembl releases
# ftp://ftp.ensemblgenomes.org/pub/release-42/metazoa/
# makeTxDbFromEnsembl(organism="Folsomia candida",
#                     release=42,
#                     circ_seqs=DEFAULT_CIRC_SEQS,
#                     server="ftp://ftp.ensemblgenomes.org/pub/release-42/metazoa/",
#                     username="anonymous", password=NULL, port=80,
#                     tx_attrib=NULL)
###########################################################################################


###########################################################################################
###########################################################################################
# Make custom TxDb using data frames
# This appears to be the best choice for nonmodel organisms; flexible and customizeable
###########################################################################################
###########################################################################################
#
# genesDF = getBM(attributes=c('ensembl_gene_id',
#                              'ensembl_transcript_id',
#                              'description',
#                              'chromosome_name',
#                              'ensembl_peptide_id'
#                               ), mart=ensembl)
#
# colnames(genesDF)[1] = 'GID'
# colnames(genesDF)[2] = 'ENSEMBLTRANS'
# colnames(genesDF)[3] = 'DESCRIPTION'
# colnames(genesDF)[4] = 'CHR'
# colnames(genesDF)[5] = 'ENSEMBLPROT'
#
# protannDF = getBM(attributes=c('ensembl_gene_id',
#                                'interpro',
#                                'pfam',
#                                'kegg_enzyme' ), mart=ensembl)
#
# colnames(protannDF)[1] = 'GID'
# colnames(protannDF)[2] = 'INTERPRO'
# colnames(protannDF)[3] = 'PFAM'
# colnames(protannDF)[4] = 'PATH'
#
# pfscanDF = getBM(attributes=c('ensembl_gene_id',
#                               'pfscan'), mart=ensembl)
#
# colnames(pfscanDF)[1] = 'GID'
# colnames(pfscanDF)[2] = 'PROSITE'
# pfscanDF_annotated_only <- pfscanDF[grep("PS", pfscanDF$PROSITE),]
#
# goDF = getBM(attributes=c('ensembl_gene_id', 'go_id', 'go_linkage_type'), mart=ensembl)
# colnames(goDF)[1] = 'GID'
# colnames(goDF)[2] = 'GO'
# colnames(goDF)[3] = 'EVIDENCE'
#
# goDF_annotated_only <- goDF[grep("GO", goDF$GO),]
#
# makeOrgPackage(gene_info=genesDF,
#                prot=protannDF,
#                prosite=pfscanDF_annotated_only,
#                go=goDF_annotated_only,
#                version="0.1",
#                author = " <>",
#                maintainer = " <>",
#                outputDir = ".",
#                tax_id="158441",
#                genus="Folsomia",
#                species="candida",
#                goTable="go")
#
# detach("package:org.Fcandida.eg.db", unload=TRUE)
# install.packages("./org.Fcandida.eg.db", repos = NULL, type="source")
# library(org.Fcandida.eg.db)
###########################################################################################


###########################################################################################
######## PRODUCE INPUT FOR BMDExpress2 ####################################################
###########################################################################################

# Create input file...
bmdexpress <- as.data.frame(lognorm.read.counts) # log2 normalized, size factor normalized, 1 added
head(bmdexpress)
bmdexpress <- cbind(SampleID=row.names(bmdexpress), bmdexpress, stringsAsFactors=F)
head(bmdexpress)
bmdexpress <- rbind( Dose=c("Dose",as.character(decoder.data$group.ID)), bmdexpress, stringsAsFactors=F)
# bmdexpress <- rbind( Dose=c("Dose","0","1","2"), bmdexpress, stringsAsFactors=F)
head(bmdexpress)
write.table(bmdexpress, file = "bmdexpress_input.txt", quote = F, sep = "\t", row.names = F, col.names = T)


# Create probe mappings... Should only need to be done once

#bmdexpressProbeMap <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id'), mart=ensembl)
# Rename columns...
#colnames(bmdexpressProbeMap)[1] = 'Array Probe'
#colnames(bmdexpressProbeMap)[2] = 'Category Component'
# Write file...
#write.table(bmdexpressProbeMap, file = "bmdexpressProbeMap.txt", quote = F, sep = "\t", row.names = F, col.names = T)
#bmdexpressCategoryMap <- getBM(attributes=c('go_id','name_1006','ensembl_transcript_id'), mart=ensembl)
# Rename columns...
#colnames(bmdexpressCategoryMap)[1] = 'Category ID'
#colnames(bmdexpressCategoryMap)[2] = 'Category Name'
#colnames(bmdexpressCategoryMap)[3] = 'Category Component'
# May need to remove empty transcript lines...
# Write file...
#write.table(bmdexpressCategoryMap, file = "bmdexpressCategoryMap.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#                       go_id                          GO term accession
#                   name_1006                               GO term name
#             definition_1006                         GO term definition
#             go_linkage_type                      GO term evidence code
#              namespace_1003                                  GO domain


###########################################################################################
###########################################################################################
#### ClusterProfiler ######################################################################
###########################################################################################
###########################################################################################

# Create geneList object from DESeq2 results
d <- data.frame(row.names(DESeqResults),DESeqResults[,2])
d2 <- d
d2[,1] <- as.character(d[,1]) # For clusterProfiler, where character is needed...
## 1st column is gene ID
## 2nd column is FC

## feature 1: numeric vector
geneList = d[,2]
## feature 2: named vector
names(geneList) = as.character(d[,1])
## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)

genesymbols=d[,1]

# columns(org.Fcandida.eg.db)
## [1] "CHR"          "DESCRIPTION"  "ENSEMBLPROT"  "ENSEMBLTRANS" "EVIDENCE"
## [6] "EVIDENCEALL"  "GID"          "GO"           "GOALL"        "INTERPRO"
## [11] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PROSITE"
# keytypes(org.Fcandida.eg.db)
## [1] "DESCRIPTION"  "ENSEMBLPROT"  "ENSEMBLTRANS" "EVIDENCE"     "EVIDENCEALL"
## [6] "GID"          "GO"           "GOALL"        "INTERPRO"     "ONTOLOGY"
## [11] "ONTOLOGYALL"  "PATH"         "PFAM"         "PROSITE"

# for (i in keytypes(org.Fcandida.eg.db)) { print(i); print(head(keys(org.Fcandida.eg.db, keytype=i))) }
# [1] "DESCRIPTION"
# [1] "(+)-caryolan-1-ol synthase"
# [2] "(+)-germacrene D synthase"
# [3] "(-)-delta-cadinene synthase"
# [4] "(3S,6E)-nerolidol synthase"
# [5] "(E)-2-epi-beta-caryophyllene synthase"
# [6] "(R)-mandelonitrile lyase 2"
# [1] "ENSEMBLPROT"
# [1] "OXA36501" "OXA36502" "OXA36503" "OXA36504" "OXA36505" "OXA36506"
# [1] "ENSEMBLTRANS"
# [1] "OXA36501" "OXA36502" "OXA36503" "OXA36504" "OXA36505" "OXA36506"
# [1] "EVIDENCE"
# [1] "IEA"
# [1] "EVIDENCEALL"
# [1] "IEA"
# [1] "GID"
# [1] "Fcan01_02664" "Fcan01_03138" "Fcan01_00050" "Fcan01_03825" "Fcan01_02931"
# [6] "Fcan01_01851"
# [1] "GO"
# [1] "GO:0000012" "GO:0000015" "GO:0000027" "GO:0000030" "GO:0000045"
# [6] "GO:0000049"
# [1] "GOALL"
# [1] "GO:0000002" "GO:0000003" "GO:0000012" "GO:0000015" "GO:0000026"
# [6] "GO:0000027"
# [1] "INTERPRO"
# [1] ""          "IPR000001" "IPR000003" "IPR000007" "IPR000008" "IPR000009"
# [1] "ONTOLOGY"
# [1] "BP" "CC" "MF"
# [1] "ONTOLOGYALL"
# [1] "BP" "CC" "MF"
# [1] "PATH"
# [1] ""               "00010+1.1.1.1"  "00010+1.1.1.27" "00010+1.8.1.4"
# [5] "00010+2.3.1.12" "00010+2.7.1.1"
# [1] "PFAM"
# [1] ""        "PF00001" "PF00002" "PF00003" "PF00004" "PF00005"
# [1] "PROSITE"
# [1] "PS01031" "PS01033" "PS01179" "PS01180" "PS01225" "PS50001"


# head(bitr(d[,1], fromType="GID", toType=c("GO"), OrgDb="org.Fcandida.eg.db"))


# MF = molecular function
# BP = biological process
# CC = cellular component
# GO is allowed as ont, but unsure what this returns

# groupGO(), enrichGO(), and gseGO() may be useful for various applications.

ggoMF_level2 <- groupGO(gene     = d2[,1],
               OrgDb    = org.Fcandida.eg.db,
               keytype  = "GID",
               ont      = "MF",
               level    = 2,
               readable = FALSE)

ggoBP_level2 <- groupGO(gene     = d2[,1],
                 OrgDb    = org.Fcandida.eg.db,
                 keytype  = "GID",
                 ont      = "BP",
                 level    = 2,
                 readable = FALSE)

ego <- gseGO(geneList     = geneList,
             keyType      = "GID",
             OrgDb        = org.Fcandida.eg.db,
             ont          = "CC",
             nPerm        = 1000,
             minGSSize    = 100,
             maxGSSize    = 500,
             pvalueCutoff = 0.05,
             verbose      = FALSE)

pdf(file="Gene_ontologies.pdf")
# shoCategory sets the number of categories to display
barplot(ggoMF_level2, drop=TRUE, showCategory=100)
barplot(ggoBP_level2, drop=TRUE, showCategory=100)
dev.off()


###########################################################################################
###########################################################################################
#### GOexpress ############################################################################
###########################################################################################
###########################################################################################
#
# allgenes.Ensembl = getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id', 'description'), mart=ensembl)
# colnames(allgenes.Ensembl)[1] = 'gene_id'
# allGO.Ensembl = getBM(attributes=c('go_id', 'name_1006', 'namespace_1003'), mart=ensembl)
# GOgenes.Ensembl = getBM(attributes=c('ensembl_gene_id', 'go_id'), mart=ensembl)
# colnames(GOgenes.Ensembl)[1] = 'gene_id'
# GOgenes.Ensembl = GOgenes.Ensembl[GOgenes.Ensembl$go_id != '',]
# GOgenes.Ensembl = GOgenes.Ensembl[GOgenes.Ensembl$gene_id != '',]
# save(GOgenes.Ensembl, file='GOgenes.Ensembl.rda')
# save(allGO.Ensembl, file='allGO.Ensembl.rda')
# save(allgenes.Ensembl, file='allgenes.Ensembl.rda')
#
# go_express_analysis <- GO_analyse(eSet=minimalSet,
#                                   f='sample.ID',
#                                   GO_genes=GOgenes.Ensembl,
#                                   all_GO=allGO.Ensembl,
#                                   all_genes=allgenes.Ensembl)
#
# BP.5 <- subset_scores(result = AlvMac_results.pVal, namespace = "biological_process",
#                       total = 5, p.val=0.05) # requires 5 or more associated genes
#
# <- subset( go_express_analysis.pVal$GO,
#            total_count >= 5 & p.val<0.05 & namespace_1003=='biological_process')
#
# heatmap_GO(go_id = "GO:0016507",
#            result = go_express_analysis,
#            eSet=minimalSet, cexRow=0.4,
#            cexCol=1,
#            cex.main=1,
#            main.Lsplit=30)
#
# expression_plot(gene_id = "Fcan01_154",
#                 result = go_express_analysis,
#                 eSet=minimalSet,
#                 x_var = "sample.ID",
#                 title.size=1.5,
#                 legend.title.size=10,
#                 legend.text.size=10,
#                 legend.key.size=15)
#
###########################################################################################

