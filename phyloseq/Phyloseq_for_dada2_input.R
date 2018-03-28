# Commands to run phyloseq visualizations and data analysis from the output of DADA2.
# Author: Matt Meier

library(phyloseq)
library(gridExtra)
library(dada2)
library(msa)
library(ggplot2)
library(plyr)
#library(phangorn) # Dependency for trees... only load during tree building
library(tidyverse)
library(DESeq2)

setwd("directory")

#seqs <- getSequences(seqtab)
#names(seqs) <- seqs
## mult <- msa(seqs, method="ClustalW", type="dna", order="input") ## Use other algorithm on server instead!
library(phangorn)
phang.align <- phyDat(as(alignment,"matrix"),type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) ## Tip order is not the same as sequence order!
fit = pml(treeNJ, data=phang.align)
# negative edges lengths changed to 0
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))

detach("package:phangorn", unload=TRUE)

sampledata.df <- read.table("./SampleData.16S.txt", sep="\t",header=T)
# row.names(seqtab.nochim) <- row.names(sampledata.df) ### Perform a sanity check on this!!!! Critical to have correct sample naming.
row.names(sampledata.df) <- gsub("\\.", "_", row.names(sampledata.df))
sampledata.df$Silver_Concentration <- as.factor(sampledata.df$Silver_Concentration)
sampleOrder=row.names(sampledata.df)

#ITS
#ps <- phyloseq(tax_table(taxa_its), sample_data(sampledata.df), otu_table(seqtab.nochim, taxa_are_rows=FALSE), phy_tree=(fitGTR$tree))
#16S
ps <- phyloseq(tax_table(taxa_silva), sample_data(sampledata.df), otu_table(seqtab.nochim, taxa_are_rows=FALSE), phy_tree=(fitGTR$tree))
ps <- prune_samples(sample_names(ps) !="GRDI-ECO_SM-LB-VAB2012-20140113-AGEZHA0-A_16S_S0054C9", ps)
table(tax_table(ps)[, "Phylum"], exclude=NULL)

ps0 <- subset_taxa(ps,  !is.na(Phylum) &! Phylum %in% c("", "uncharacterized"))
table(tax_table(ps0)[, "Phylum"], exclude=NULL)

prevdf = apply(X = otu_table(ps0), MARGIN = ifelse(taxa_are_rows(ps0), yes =1, no =2), FUN = function(x) {sum(x>0)})
prevdf = data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(ps0),tax_table(ps0))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
filterPhyla = c("Woesearchaeota", "Microgenomates", "Euryarchaeota", "Diapherotrites")
ps1=subset_taxa(ps0, !Phylum %in% filterPhyla)


# Filter by prevalence

prevdf1 = subset(prevdf, Phylum  %in% get_taxa_unique(ps1,"Phylum"))
ggplot(prevdf1,  aes(TotalAbundance, Prevalence/nsamples(ps1), color=Phylum)) + 
  geom_hline(yintercept=0.01, alpha=0.5, linetype=2) +
  geom_point(size=2, alpha=0.7) + scale_x_log10()+  xlab("Total Abundance") + 
  ylab("Prevalence [Frac. Samples]")+ facet_wrap(~ Phylum) + theme(legend.position="none")

prevalenceThreshold = 0.01*nsamples(ps1)
keepTaxa=rownames(prevdf1)[(prevdf1$Prevalence>=prevalenceThreshold)]
ps2=prune_taxa(keepTaxa, ps1)
length(get_taxa_unique(ps2,taxonomic.rank="Genus"))

### Plot Richness

plot_richness(ps0, x="Silver_Concentration", measures=c("Observed", "Shannon", "Simpson", "Chao1"), color="Silver_Concentration") + theme_bw() + geom_point(size=3) +  labs(x = "Silver Concentration", colour="Silver Concentration") + geom_boxplot()

### Plot trees

ps3=tax_glom(ps2,"Genus" , NArm=TRUE)
h1=0.2
ps4=tip_glom(ps2, h = h1)

multiPlotTitleTextSize=12
p2tree= plot_tree(ps2,method="treeonly",
                  ladderize="left",
                  title="Before Agglomeration") + 
  theme(plot.title=element_text(size= multiPlotTitleTextSize))
p3tree= plot_tree(ps3,method="treeonly",
                  ladderize="left",
                  title="By Genus") +
  theme(plot.title=element_text(size= multiPlotTitleTextSize))
p4tree=plot_tree(ps4,method="treeonly",
                 ladderize="left",
                 title="By Height") +
  theme(plot.title=element_text(size=multiPlotTitleTextSize))
# group plots together
grid.arrange(nrow=1, p2tree, p3tree, p4tree)


plot_abundance=function(physeq,title="",Facet="Order",Color="Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f= subset_taxa(physeq, Phylum %in% c("p__Ascomycota"))
  mphyseq=psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance> 0)
  ggplot(data=mphyseq,mapping=aes_string(x="Silver_Concentration",y="Abundance",color= Color,fill= Color)) +
    geom_boxplot(fill=NA)+geom_point(size=1,alpha=0.3,position=position_jitter(width=0.3)) +
    facet_wrap(facets= Facet) + scale_y_log10()+theme(legend.position="none")
}

ps3ra=transform_sample_counts(ps3, function(x){x/sum(x)})
ps2ra=transform_sample_counts(ps2, function(x){x/sum(x)})

plotBefore=plot_abundance(ps3,"")
plotAfter=plot_abundance(ps3ra,"")

grid.arrange(nrow=2, plotBefore, plotAfter)

plot_abundance(ps1,Facet="Genus",Color=NULL)
plot_abundance(ps2,Facet="Genus",Color=NULL)
plot_abundance(ps2ra,Facet="Genus",Color=NULL)

### ORDINATION
ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
plot_ordination(ps, ord.nmds.bray, color="Silver_Concentration", title="Bray NMDS") + geom_point(size=5) + labs(colour="Silver Concentration")
ord.nmds.jaccard <- ordinate(ps, method="NMDS", distance="jaccard")
plot_ordination(ps, ord.nmds.jaccard, color="Silver_Concentration", title="Jaccard NMDS") + geom_point(size=5) + labs(colour="Silver Concentration")
ord.nmds.wunifrac <- ordinate(ps, method="NMDS", distance="wunifrac")
plot_ordination(ps, ord.nmds.wunifrac, color="Silver_Concentration", title="Weighted Unifrac NMDS")
#ord.rda.wunifrac <- ordinate(ps, method="RDA", distance="wunifrac")
#plot_ordination(ps, ord.rda.wunifrac, color="Silver_Concentration", title="Weighted Unifrac RDA")
#ord.rda.unifrac <- ordinate(ps, method="RDA", distance="unifrac")
#plot_ordination(ps, ord.rda.unifrac, color="Silver_Concentration", title="Unifrac RDA")
ord.dca.wunifrac <- ordinate(ps, method="DCA", distance="wunifrac")
plot_ordination(ps, ord.dca.wunifrac, color="Silver_Concentration", title="Weighted Unifrac DCA") + geom_point(size=5) + labs(colour="Silver Concentration")
ord.dca.bray <- ordinate(ps, method="DCA", distance="bray")
plot_ordination(ps, ord.dca.bray, color="Silver_Concentration", title="Bray DCA")

ps_glom_species <- tax_glom(ps, taxrank="Species", NArm=FALSE)
ps_glom_species_transformed <-transform_sample_counts(ps_glom_species, function(x) x / sum(x))

plot_tree(ps_glom_species_transformed, method="sampledodge", size="Abundance", justify="yes please", ladderize="left" , color="Silver_Concentration") +
       scale_size_continuous(range = c(1, 3))

library(RColorBrewer)

plot_bar(ps_glom_species_transformed, fill="Genus")

filterfun1 = function(x){
  xprop = (x / sum(x))
  xprop[xprop < (1e-5)] <- 0
  return(xprop)
}

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- transform_sample_counts(ps, function(x) 100*(x/sum(x)))
ps.top20 <- prune_taxa(top20, ps.top20)

getPalette = colorRampPalette(brewer.pal(11, "Spectral"))

getPalette = colorRampPalette(brewer.pal(9, "Set1")) # not bad - but may need to rearrange
getPalette = colorRampPalette(brewer.pal(8, "Set2")) # too muted
getPalette = colorRampPalette(brewer.pal(8, "Set3"))
getPalette = colorRampPalette(brewer.pal(8, "Accent")) # Ugly
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

plot_bar(ps.top20, x="Total_treatment", fill="Genus") + 
  labs(x="Silver Concentration") + 
  facet_grid(~Silver_Concentration,scales="free") +
  scale_fill_manual(values=getPalette(22))

top100 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:100]
ps.top100 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top100 <- prune_taxa(top100, ps.top100)
plot_bar(ps.top100, x="Sample", fill="Genus") + labs(x="Silver Concentration") + facet_grid(~Silver_Concentration,scales="free")

### Kostic example, adapted for GRDI data from dada2
kostic <- ps
kostic <- prune_samples(sample_sums(kostic) > 500, kostic)
diagdds = phyloseq_to_deseq2(kostic, ~ Silver_Concentration)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
res = results(diagdds)

res60 <- results(diagdds, contrast=c("Silver_Concentration","60","0"), cooksCutoff = FALSE)
res145 <- results(diagdds, contrast=c("Silver_Concentration","145","0"), cooksCutoff = FALSE)
res347 <- results(diagdds, contrast=c("Silver_Concentration","347","0"), cooksCutoff = FALSE)
res833 <- results(diagdds, contrast=c("Silver_Concentration","833","0"), cooksCutoff = FALSE)
res2000 <- results(diagdds, contrast=c("Silver_Concentration","2000","0"), cooksCutoff = FALSE)

alpha = 0.01

sigtab60 = res60[which(res60$padj < alpha), ]
sigtab145 = res145[which(res145$padj < alpha), ]
sigtab347 = res347[which(res347$padj < alpha), ]
sigtab833 = res833[which(res833$padj < alpha), ]
sigtab2000 = res2000[which(res2000$padj < alpha), ]

sigtab60 = cbind(as(sigtab60, "data.frame"), as(tax_table(ps)[rownames(sigtab60), ], "matrix"))
sigtab145 = cbind(as(sigtab145, "data.frame"), as(tax_table(ps)[rownames(sigtab145), ], "matrix"))
sigtab347 = cbind(as(sigtab347, "data.frame"), as(tax_table(ps)[rownames(sigtab347), ], "matrix"))
sigtab833 = cbind(as(sigtab833, "data.frame"), as(tax_table(ps)[rownames(sigtab833), ], "matrix"))
sigtab2000 = cbind(as(sigtab2000, "data.frame"), as(tax_table(ps)[rownames(sigtab2000), ], "matrix"))


# Plot results
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

sigtableList <- NULL
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

combined_sigtableListNA.RM <- do.call(rbind, sigtableListNA.RM)

#ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
#  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

## Genus plots, no N/A
ggplot(subset(sigtab833, !is.na(Genus)), aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Make phyloseq object from the differentially abundant genera
dataSamplesPruned <- prune_samples(sample_sums(ps)>=1000, ps)
dataSamplesPrunedRarefied <- rarefy_even_depth(dataSamplesPruned)
dataSamplesPrunedRarefiedSpeciesGlom <- tax_glom(dataSamplesPrunedRarefied, taxrank="Species", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
dataSamplesPrunedRarefiedSpeciesGlom <- tax_glom(dataSamplesPrunedRarefied, taxrank="Species", NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))
dataSamplesPrunedRarefiedGenusGlom <- tax_glom(dataSamplesPrunedRarefied, taxrank="Genus", NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))

plot_heatmap(dataSamplesPrunedRarefiedSpeciesGlom, "RDA", "unifrac", "Silver_Concentration", "Genus",    sample.order=row.names(sampledata.df))

dataSamplesSpeciesGlomRarefied <- tax_glom(dataSamplesPrunedRarefied, taxrank="Species", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
# genera_affected <- subset_taxa(ps, Genus %in% levels(subset(sigtab, !is.na(Genus))$Genus))
genera_affected_sp_glom <- subset_taxa(dataSamplesPrunedRarefiedSpeciesGlom, Genus %in% levels(subset(combined_sigtableListNA.RM, !is.na(Genus))$Genus))
genera_affected_sp_glom_rarefied <- subset_taxa(dataSamplesSpeciesGlomRarefied, Genus %in% levels(subset(combined_sigtableListNA.RM, !is.na(Genus))$Genus))

plot_heatmap(genera_affected_sp_glom, "NMDS", "jaccard", "Silver_Concentration", "Genus", sample.order=sampleOrder) ## Nice

plot_heatmap(genera_affected_sp_glom, "RDA", "none", "Silver_Concentration", "Genus", sample.order=sampleOrder)
plot_heatmap(genera_affected_sp_glom_rarefied, "MDS", "unifrac", "Silver_Concentration", "Genus",    sample.order=sampleOrder)
plot_heatmap(dataSamplesPrunedRarefiedSpeciesGlom, "RDA", "unifrac", "Silver_Concentration", "Genus",    sample.order=sampleOrder)

# ORDERING TAXA
sigtab2000$taxID <- row.names(sigtab2000)
sigtab833$taxID <- row.names(sigtab833)
sigtab347$taxID <- row.names(sigtab347)
sigtab145$taxID <- row.names(sigtab145)
sigtab60$taxID <- row.names(sigtab60)
combined_sigtableListNA.RM$taxID <- row.names(combined_sigtableListNA.RM)

orderedGenera833 <- as.character(arrange(sigtab833, log2FoldChange)$taxID[!is.na(arrange(sigtab833, log2FoldChange)$taxID)])
orderedGenera60 <- as.character(arrange(sigtab60, log2FoldChange)$taxID[!is.na(arrange(sigtab60, log2FoldChange)$taxID)])
orderedGenera145 <- as.character(arrange(sigtab145, log2FoldChange)$taxID[!is.na(arrange(sigtab145, log2FoldChange)$taxID)])
orderedGenera347 <- as.character(arrange(sigtab347, log2FoldChange)$taxID[!is.na(arrange(sigtab347, log2FoldChange)$taxID)])
orderedGenera2000 <- as.character(arrange(sigtab2000, log2FoldChange)$taxID[!is.na(arrange(sigtab2000, log2FoldChange)$taxID)])

orderedGeneraAll <- as.character(arrange(combined_sigtableListNA.RM, log2FoldChange)$taxID[!is.na(arrange(combined_sigtableListNA.RM, log2FoldChange)$taxID)])
orderedGeneraByPval <- as.character(arrange(combined_sigtableListNA.RM, padj)$taxID[!is.na(arrange(combined_sigtableListNA.RM, log2FoldChange)$taxID)])

### ADD THIS IN PLACE FOR THE HEATMAP TO ELIMINATE UNASSIGNED GENERA
# subset_taxa(dataSamplesPrunedRarefied,  !is.na(Genus)

plot_heatmap(subset_taxa(dataSamplesPrunedRarefied,  !is.na(Genus)), "none", "none", "Silver_Concentration", 
                  taxa.label="Genus", sample.order=sampleOrder, taxa.order=orderedGeneraAll, max.label=5000) +
                  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10)) +
                  labs(x="Silver Concentration")

nitrifying_bacteria <- c("Nitrosomonas", "Nitrosococcus", "Nitrobacter", "Nitrococcus", "Nitrospina", "Nitrospira", "Nitrosophaera")
ps_nitrifying <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
#ps_nitrifying <- subset_taxa(ps, Genus %in% nitrifying_bacteria)
ps_nitrifying <- subset_taxa(ps_nitrifying, Genus %in% nitrifying_bacteria)
ps_nitrifying_glom <- tax_glom(ps_nitrifying,taxrank="Genus")

plot_heatmap(ps_nitrifying_glom,  "Silver_Concentration", 
             taxa.label="Genus", sample.order=sampleOrder,  max.label=5000) +
             theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10)) +
             labs(x="Silver Concentration")
plot_bar(ps_nitrifying_glom, fill="Genus")
plot_bar(ps_nitrifying_glom, x="Silver_Concentration", fill="Genus")

taxdf <- as.data.frame(tax_table(ps))
taxdf$taxid <- row.names(taxdf)
ps_transformed <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
nitro_tax <- subset(taxdf, grepl("Nitro", taxdf$Kingdom) | grepl("Nitro", taxdf$Phylum) | grepl("Nitro", taxdf$Class) | grepl("Nitro", taxdf$Order) | grepl("Nitro", taxdf$Genus) | grepl("Nitro", taxdf$Species) )
nitro_tax_ids <- row.names(nitro_tax)
ps_nitro=prune_taxa(nitro_tax_ids, ps_transformed)
ps_nitro_glom <- tax_glom(ps_nitro, taxrank="Genus")
plot_bar(ps_nitro_glom, x="Silver_Concentration", fill="Genus")
plot_heatmap(ps_nitro,  "Silver_Concentration", 
             taxa.label="Genus", sample.order=sampleOrder,  max.label=5000) +
             theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10)) +
             labs(x="Silver Concentration")

