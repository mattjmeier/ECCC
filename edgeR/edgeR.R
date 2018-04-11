
# Run edgeR on a Biom file

# First make sure you have done the following in Qiime:

# source /qiime1.9/bin/activate
# summarize_taxa.py -i ./otu_table.biom -L 6 -a -o taxa-absolute # GENUS LEVEL... L7 for species
# biom  convert -i otu_table_L6.biom -o otu_table_L6_json.biom --to-json

# install.packages(c('biom','vegan'), repo='http://cran.wustl.edu')


source("https://bioconductor.org/biocLite.R")
biocLite("biomformat")
library(biomformat)
biocLite("edgeR")
library("edgeR")

source(wrapper) ## This is the edgeR wrapper written by Dan Knights

setwd("directory")

genus.biom <- read_biom("./otu_table_L6_json.biom")
species.biom <- read_biom("./otu_table_L7_json.biom")

genus <- as.matrix(biom_data(genus.biom))
genus <- t(genus)

species <- as.matrix(biom_data(species.biom))
species <- t(species)

map <- read.table('SampleData.txt', sep='\t', head=T, row.names=1)

common.ids <- intersect(rownames(map), rownames(genus))
common.ids.sp <- intersect(rownames(map), rownames(species))

genus <- genus[common.ids,]
species <- species[common.ids,]
map <- map[common.ids,]

dim(genus)

genus <- genus[,colMeans(genus > 0) >= 0.1]
species <- species[,colMeans(species > 0) >= 0.1]

dim(genus)

colnames(genus)[1:10]
# colnames(genus) <- sapply(strsplit(colnames(genus),'p__'),'[',2)
# colnames(species) <- sapply(strsplit(colnames(species),'p__'),'[',2)

colnames(genus) <- sapply(strsplit(colnames(genus),';'),function(xx) paste(paste(substr(xx[-c(1,length(xx))],4,7),collapse=';'),substring(xx[length(xx)],4),sep=';'))
                                                                                               
results.genus <- glm.edgeR(x=map$Treatment, Y=genus, covariates=map$Replicate) ## , covariates=map.16S$covariate

topTags(results.genus)

QL.results.genus.group.replicate <- glmQLF.edgeR(x=map$group, Y=genus,  covariates=map$Replicate)
QL.results.species.group.replicate <- glmQLF.edgeR(x=map$group, Y=species,  covariates=map$Replicate)

write.table(topTags(results, n=Inf)$table, file='edgeR_results.txt', sep="\t",quote=F, col.names=NA)
