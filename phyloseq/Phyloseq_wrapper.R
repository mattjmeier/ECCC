# Phyloseq and differential abundance analysis starting with a biom file 

library(phyloseq)
library(ggplot2)
library(plyr)
library(tidyverse)
library("DESeq2")


## filter(as.data.frame(tax_table(data)), str_detect(Genus, "Rhodanobacter"))


setwd("directory")
# files <- list.files(path="./", pattern="*.biom$", full.names=T, recursive=FALSE)

#bottomRank<-paste("Rank6")
#topRank<-paste("Rank2")

x = read_tree_greengenes("97_otus.tree")
data<-import_biom("./otu_table.biom", parseFunction=parse_taxonomy_greengenes)
sampledata<-sample_data(read.table("SampleData.txt", sep="\t",header=T))
sampledata$Silver.Concentration <- as.factor(sampledata$Silver.Concentration)
sampledata$Total_treatment <- factor(sampledata$Total_treatment, levels = c("None 0", "HA 0", "PVP 0", "Silver 60", "Silver 145", "Silver 347", "Silver 833", "Silver 2000"))
dataSamples<-merge_phyloseq(data,sampledata)
dataSamples<-merge_phyloseq(dataSamples, x)

### Phylogenetic Trees ###
mdt = fast_melt(dataSamples)
prevdt = mdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]
keepTaxa = prevdt[(Prevalence >= 10 & TotalCounts > 3), TaxaID]
dataKeptTaxa = prune_taxa(keepTaxa, dataSamples)
tipg = tip_glom(dataKeptTaxa, h = 0.05)
taxg = tax_glom(dataKeptTaxa, taxrank="Genus", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

# Transform to relative abundance
#tipg <- transform_sample_counts(tipg, function(x) x / sum(x))
taxg <- transform_sample_counts(taxg, function(x) x / sum(x))

#ntaxa(tipg)
ntaxa(taxg)

# Plot tree
plot_tree(taxg, size = "Abundance", method = "sampledodge",
          color = "Silver.Concentration",
          justify = "yes please", 
          ladderize = "left") +
          scale_size_continuous(range = c(1, 3))

dist_methods <- unlist(distanceMethodList)
print(dist_methods)
dist_methods = dist_methods[-which(dist_methods=="ANY")]
dist_methods = dist_methods[-which(dist_methods=="manhattan")]
print(dist_methods)
plist <- vector("list", length(dist_methods)) #### LENGTH MUST AGREE WITH WHAT YOU DO BELOW

for( i in dist_methods ){
  iDist <- distance(dataSamples, method=i)
  iMDS  <- ordinate(dataSamples, "MDS", distance=iDist)
  p <- NULL
  p <- plot_ordination(dataSamples, iMDS, color="Silver.Concentration", shape="Amendment")
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  plist[[i]] = p
}

df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=Treatment, shape=Body_site))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for AgNP / Antibiotics Microbiome Dataset")
p


df_sub<-df[df$distance %in% c("-1","bray","canberra","cao","cc","co","dpcoa","g","hk","jaccard","jsd","kulcynski","l","t","unifrac","z"),]
names(df_sub)[1] <- "distance"
p = ggplot(df_sub, aes(Axis.1, Axis.2, color=Antibiotics, shape=Silver_Treatment))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for AgNP / Antibiotics Microbiome Dataset")
p

### DESeq2 Differential Counts

dataSamplesDESeq2 <- phyloseq_to_deseq2(dataSamples, ~ Silver_Concentration)
dataSamplesDESeq2_result <- DESeq(dataSamplesDESeq2, test="Wald", fitType="parametric") ### This probably won't work for denoised data, e.g., dada2 results
res <- results(dataSamplesDESeq2_result, cooksCutoff = FALSE)
res60 <- results(dataSamplesDESeq2_result, contrast=c("Silver.Concentration","60","0"), cooksCutoff = FALSE) ## Specify what to compare!
res145 <- results(dataSamplesDESeq2_result, contrast=c("Silver.Concentration","145","0"), cooksCutoff = FALSE)
res347 <- results(dataSamplesDESeq2_result, contrast=c("Silver.Concentration","347","0"), cooksCutoff = FALSE)
res833 <- results(dataSamplesDESeq2_result, contrast=c("Silver.Concentration","833","0"), cooksCutoff = FALSE)
res2000 <- results(dataSamplesDESeq2_result, contrast=c("Silver.Concentration","2000","0"), cooksCutoff = FALSE)

alpha = 0.01
sigtab = res[which(res$padj < alpha), ]

sigtab60 = res60[which(res60$padj < alpha), ]
sigtab145 = res145[which(res145$padj < alpha), ]
sigtab347 = res347[which(res347$padj < alpha), ]
sigtab833 = res833[which(res833$padj < alpha), ]
sigtab2000 = res2000[which(res2000$padj < alpha), ]

sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dataSamples)[rownames(sigtab), ], "matrix"))

sigtab60 = cbind(as(sigtab60, "data.frame"), as(tax_table(dataSamples)[rownames(sigtab60), ], "matrix"))
sigtab145 = cbind(as(sigtab145, "data.frame"), as(tax_table(dataSamples)[rownames(sigtab145), ], "matrix"))
sigtab347 = cbind(as(sigtab347, "data.frame"), as(tax_table(dataSamples)[rownames(sigtab347), ], "matrix"))
sigtab833 = cbind(as(sigtab833, "data.frame"), as(tax_table(dataSamples)[rownames(sigtab833), ], "matrix"))
sigtab2000 = cbind(as(sigtab2000, "data.frame"), as(tax_table(dataSamples)[rownames(sigtab2000), ], "matrix"))

head(sigtab145)

# Plot results
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

sigtableList <- list(sigtab60,sigtab145,sigtab347,sigtab833,sigtab2000)

for (sigtable in sigtableList) {

# Phylum order
x = tapply(sigtable$log2FoldChange, sigtable$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtable$Phylum = factor(as.character(sigtable$Phylum), levels=names(x))
# Genus order
x = tapply(sigtable$log2FoldChange, sigtable$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtable$Genus = factor(as.character(sigtable$Genus), levels=names(x))
 
}
  
sigtableListNA.RM <-lapply(sigtableList, subset, !is.na(Genus))
for (sigtable in sigtableListNA.RM) {
  
  sigtable$Genus
  
}

combined_sigtableListNA.RM <- do.call(rbind, sigtableListNA.RM)


#ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
#  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

## Genus plots, no N/A
ggplot(subset(sigtab, !is.na(Genus)), aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
   theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Make phyloseq object from the differentially abundant genera
dataSamplesPruned <- prune_samples(sample_sums(dataSamples)>=1000, dataSamples)
dataSamplesPrunedRarefied <- rarefy_even_depth(dataSamplesPruned)
dataSamplesPrunedRarefiedSpeciesGlom <- tax_glom(dataSamplesPrunedRarefied, taxrank="Species", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
dataSamplesPrunedRarefiedSpeciesGlom <- tax_glom(dataSamplesPrunedRarefied, taxrank="Species", NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))

plot_heatmap(dataSamplesPrunedRarefiedSpeciesGlom, "RDA", "unifrac", "Silver.Concentration", "Genus",    sample.order=sampleOrder)

#dataSamplesSpeciesGlomRarefied <- tax_glom(dataSamplesRarefied, taxrank="Species", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
# genera_affected <- subset_taxa(dataSamples, Genus %in% levels(subset(sigtab, !is.na(Genus))$Genus))
genera_affected_sp_glom <- subset_taxa(dataSamplesPrunedRarefiedSpeciesGlom, Genus %in% levels(subset(combined_sigtableListNA.RM, !is.na(Genus))$Genus))
genera_affected_sp_glom_rarefied <- subset_taxa(dataSamplesSpeciesGlomRarefied, Genus %in% levels(subset(combined_sigtableListNA.RM, !is.na(Genus))$Genus))

plot_heatmap(genera_affected_sp_glom, "NMDS", "jaccard", "Silver.Concentration", "Genus", sample.order=sampleOrder) ## Nice

plot_heatmap(genera_affected_sp_glom, "RDA", "none", "Silver.Concentration", "Genus", sample.order=sampleOrder)
plot_heatmap(genera_affected_sp_glom_rareified, "MDS", "unifrac", "Silver.Concentration", "Genus", sample.order=sampleOrder)
plot_heatmap(dataSamplesPrunedRarefiedSpeciesGlom, "RDA", "unifrac", "Silver.Concentration", "Genus", sample.order=sampleOrder)


# ORDERING TAXA
sigtab833$taxID <- row.names(sigtab833)
sigtab60$taxID <- row.names(sigtab60)
orderedGenera <- as.character(arrange(sigtab833, log2FoldChange)$taxID[!is.na(arrange(sigtab833, log2FoldChange)$taxID)])
orderedGenera60 <- as.character(arrange(sigtab60, log2FoldChange)$taxID[!is.na(arrange(sigtab60, log2FoldChange)$taxID)])

orderedGeneraFiltered <- intersect(taxa_names(genera_affected_sp_glom), orderedGenera)

plot_heatmap(genera_affected_sp_glom, "MDS", "jaccard", "Silver.Concentration", "Genus",  sample.order=sampleOrder, taxa.order= orderedGenera  )
plot_heatmap(dataSamples, "none", "none", "Silver.Concentration", taxa.label="Genus", sample.order=sampleOrder, taxa.order= orderedGenera, max.label=2500)

plot_heatmap(dataSamplesPruned, "none", "none", "Silver.Concentration", taxa.label="Genus", sample.order=sampleOrder, taxa.order= orderedGenera, max.label=2500)

p <- plot_heatmap(dataSamplesPrunedRarefied, "none", "none", "Silver.Concentration", taxa.label="Genus", sample.order=sampleOrder, taxa.order= orderedGenera60, max.label=2500)
p$scales$scales[[1]]$name <- "Silver Concentration"
p + theme (axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))

data_rarefied<-rarefy_even_depth(data, sample.size = 100)
plot_bar(data_rarefied, "Phylum", fill="Genus", facet_grid=~Sample) + geom_bar( stat="identity", position="stack")


for (currentFile in files) {
  
data<-import_biom(currentFile, parseFunction=parse_taxonomy_greengenes)

colnames(tax_table(data)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

data.ord<-ordinate(data)
data.ord.pcoa <-ordinate(data, method="PCoA")
topsp <- names(sort(taxa_sums(data), TRUE)[1:300])
data.topsp <- prune_taxa(topsp, data) 

pdf(file=paste(currentFile,".", bottomRank, ".pdf", sep=""), useDingbats=FALSE)

print(plot_bar(data.topsp, topRank, fill=bottomRank, facet_grid=~Sample) 
      + geom_bar( stat="identity", position="stack")
      + theme(legend.position="bottom", legend.text=element_text(colour="blue", size=5))
      + theme(legend.key.size = unit(0.5,"line"))
      )

#### Renuka's Project

v3="./Auto_user_SN2-12-Renuka_STAGE__Cleaning_Products_16S_Metagenomic_63_037.silva.nr_v128.v3.phylip.opti_mcc.0.03.norm.0.03.biom"
v6="./Auto_user_SN2-12-Renuka_STAGE__Cleaning_Products_16S_Metagenomic_63_037.silva.nr_v128.v6.phylip.opti_mcc.0.03.norm.0.03.biom"
v3data<-import_biom(v3, parseFunction=parse_taxonomy_greengenes)
v6data<-import_biom(v6, parseFunction=parse_taxonomy_greengenes)

colnames(v3data@otu_table@.Data)
sampleOrder <- readClipboard()

print(plot_heatmap(physeq=data.topsp, method = "Sample"))
# plot_heatmap(physeq=data.topsp, "RDA", "none", "Silver.Concentration", "Genus")
# plot_heatmap(physeq=rarefied.data.topsp, "RDA", "none", "Silver.Concentration", "Family")
# plot_heatmap(physeq=data.topsp, method="RDA", distance="none", sample.label="Amendment", taxa.label="Order", sample.order = sampleOrder, taxa.order="Genus")
# plot_heatmap(data.topsp.rarefied, "NMDS", "jaccard", "Amendment", "Family", sample.order=sampleOrder)

# plot_heatmap(genera_affected_ITS.rarefied,  "NMDS", "jaccard", "group", "Family", sample.order= sampleOrder) # To produce taxa names, the number of taxa IN THE PHYLOSEQ OBJECT must be less than the max allowed (can be specified)

# plot_heatmap(genera_affected_ITS.rarefied,  "RDA", "none", "group", "Family", sample.order= sampleOrder, taxa.order="Family")

print(plot_ordination(data, data.ord))
# plot_ordination(dataSamples, ordinate(dataSamples, "CCA"), shape="Amendment", color="Silver.Concentration")
# + geom_polygon(aes(fill = "Silver.Concentration"))

dev.off()

}

## Fast melt function - possibly now built into phyloseq as psmelt(), see https://github.com/joey711/phyloseq/issues/418

fast_melt = function(physeq,
                     includeSampleVars = character(),
                     omitZero = FALSE){
  require("phyloseq")
  require("data.table")
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "TaxaID")
  # Enforce character TaxaID key
  otudt[, TaxaIDchar := as.character(TaxaID)]
  otudt[, TaxaID := NULL]
  setnames(otudt, "TaxaIDchar", "TaxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "TaxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  if(omitZero){
    # Omit zeroes and negative numbers
    mdt <- mdt[count > 0]
  }
  # Omit NAs
  mdt <- mdt[!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "TaxaID")
    # Enforce character TaxaID key
    taxdt[, TaxaIDchar := as.character(TaxaID)]
    taxdt[, TaxaID := NULL]
    setnames(taxdt, "TaxaIDchar", "TaxaID")
    # Join with tax table
    setkey(taxdt, "TaxaID")
    setkey(mdt, "TaxaID")
    mdt <- taxdt[mdt]
  }
  # includeSampleVars = c("DaysSinceExperimentStart", "SampleType")
  # includeSampleVars = character()
  # includeSampleVars = c()
  # includeSampleVars = c("aksjdflkas") 
  wh.svars = which(sample_variables(physeq) %in% includeSampleVars)
  if( length(wh.svars) > 0 ){
    # Only attempt to include sample variables if there is at least one present in object
    sdf = as(sample_data(physeq), "data.frame")[, wh.svars, drop = FALSE]
    sdt = data.table(sdf, keep.rownames = TRUE)
    setnames(sdt, "rn", "SampleID")
    # Join with long table
    setkey(sdt, "SampleID")
    setkey(mdt, "SampleID")
    mdt <- sdt[mdt]
  }
  setkey(mdt, "TaxaID")
  return(mdt)
}

# Sample names for convenience

GRDI-ECO.SM-LB-VAB2012-20130713-AGBRHA0-A.ITS.S0054CC
GRDI-ECO.SM-LB-VAB2012-20130713-AGBRHA0-B.ITS.S0054D8
GRDI-ECO.SM-LB-VAB2012-20130713-AGBRHA0-C.ITS.S005485
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZHA0-A.ITS.S0054CF
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZHA0-C.ITS.S005488
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZHA0-B.ITS.S0054DB
GRDI-ECO.SM-LB-VAB2012-20130713-AGBRCN0-A.ITS.S005484
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZCN0-A.ITS.S005487
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZCN0-C.ITS.S00549F
GRDI-ECO.SM-LB-VAB2012-20130713-AGBRCN0-C.ITS.S00549C
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZCN0-B.ITS.S005493
GRDI-ECO.SM-LB-VAB2012-20130713-AGBRCN0-B.ITS.S005490
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZPV0-B.ITS.S0054B7
GRDI-ECO.SM-LB-VAB2012-20130713-AGBRPV0-C.ITS.S0054C0
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZPV0-A.ITS.S0054AB
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZPV0-C.ITS.S0054C3
GRDI-ECO.SM-LB-VAB2012-20130713-AGBRPV0-A.ITS.S0054A8
GRDI-ECO.SM-LB-VAB2012-20130713-AGBRPV0-B.ITS.S0054B4
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ60-C.ITS.S0054AC
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ60-A.ITS.S005494
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ60-B.ITS.S0054A0
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR60-C.ITS.S0054A9
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR60-A.ITS.S005491
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR60-B.ITS.S00549D
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ145-B.ITS.S0054C4
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR145-A.ITS.S0054B5
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR145-B.ITS.S0054C1
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ145-C.ITS.S0054D0
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR145-C.ITS.S0054CD
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ145-A.ITS.S0054B8
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR347-C.ITS.S005492
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR347-A.ITS.S0054D9
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR347-B.ITS.S005486
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ347-A.ITS.S0054DC
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ347-B.ITS.S005489
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ347-C.ITS.S005495
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR833-A.ITS.S00549E
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR833-B.ITS.S0054AA
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR833-C.ITS.S0054B6
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ833-C.ITS.S0054B9
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ833-B.ITS.S0054AD
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ833-A.ITS.S0054A1
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR2000-B.ITS.S0054CE
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR2000-C.ITS.S0054DA
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ2000-A.ITS.S0054C5
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR2000-A.ITS.S0054C2
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ2000-B.ITS.S0054D1
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ2000-C.ITS.S0054DD


GRDI-ECO.SM-LB-VAB2012-20130713-AGBRCN0-B.16S.S00548A
GRDI-ECO.SM-LB-VAB2012-20130713-AGBRCN0-A.16S.S00547E
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZPV0-C.16S.S0054BD
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZHA0-B.16S.S0054D5
GRDI-ECO.SM-LB-VAB2012-20130713-AGBRHA0-C.16S.S00547F
GRDI-ECO.SM-LB-VAB2012-20130713-AGBRHA0-B.16S.S0054D2
GRDI-ECO.SM-LB-VAB2012-20130713-AGBRPV0-B.16S.S0054AE
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZHA0-A.16S.S0054C9
GRDI-ECO.SM-LB-VAB2012-20130713-AGBRPV0-C.16S.S0054BA
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZPV0-A.16S.S0054A5
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZPV0-B.16S.S0054B1
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZHA0-C.16S.S005482
GRDI-ECO.SM-LB-VAB2012-20130713-AGBRCN0-C.16S.S005496
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZCN0-C.16S.S005499
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZCN0-A.16S.S005481
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZCN0-B.16S.S00548D
GRDI-ECO.SM-LB-VAB2012-20130713-AGBRPV0-A.16S.S0054A2
GRDI-ECO.SM-LB-VAB2012-20130713-AGBRHA0-A.16S.S0054C6
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ60-C.16S.S0054A6
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ60-A.16S.S00548E
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR60-C.16S.S0054A3
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR60-B.16S.S005497
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR60-A.16S.S00548B
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ60-B.16S.S00549A
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ145-B.16S.S0054BE
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR145-A.16S.S0054AF
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR145-B.16S.S0054BB
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ145-A.16S.S0054B2
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR145-C.16S.S0054C7
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ145-C.16S.S0054CA
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR347-B.16S.S005480
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ347-C.16S.S00548F
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ347-B.16S.S005483
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR347-A.16S.S0054D3
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ347-A.16S.S0054D6
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR347-C.16S.S00548C
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR833-C.16S.S0054B0
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ833-B.16S.S0054A7
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ833-C.16S.S0054B3
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR833-A.16S.S005498
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR833-B.16S.S0054A4
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ833-A.16S.S00549B
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR2000-B.16S.S0054C8
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ2000-C.16S.S0054D7
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ2000-B.16S.S0054CB
GRDI-ECO.SM-LB-VAB2012-20140113-AGEZ2000-A.16S.S0054BF
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR2000-C.16S.S0054D4
GRDI-ECO.SM-LB-VAB2012-20130713-AGBR2000-A.16S.S0054BC



#find percent abundance of phyloseq object
obj2 <- transform_sample_counts(obj1, function(x) 100 * x/sum(x))
#melt the phyloseq object with psmelt
obj3 <- psmelt(obj2)
#rename NA's in Class level to "Unknown"
obj3$Class[is.na(obj3$Class)]="Unknown"
#make data a data.frame
obj4 <- data.frame(obj3)

class <- ggplot(obj4, aes(x=Sample, y=Abundance, fill=Class, order=as.factor(Class)))
class <- class + geom_bar()
class <- class + theme()
class

