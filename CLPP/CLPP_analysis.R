## To read in parts of excel sheet for CLPP data and run PCoA on datasets

library("xlsx")
library("ggplot2")
library("vegan")
library("phyloseq")
library("factoextra")

excelFiles <- c("./Feb 12.xls","./Feb 26.xls", "./March 12.xls", "./March 19.xls")
#excelFiles <- c("./March 19.xls")

compoundNames<-as.vector(read.table("./compound_names.txt", sep="\n"))

myDataSummary <- NULL

for (timepoint in 1:length(excelFiles)) {

workbook <- excelFiles[timepoint]

number_of_sheets=7

for (idx in 1:number_of_sheets) {

print(idx)

  
  # Raw Data
  
   # R1 <- as.vector(t(read.xlsx(workbook, sheetIndex = idx, rowIndex = 68:75, colIndex = 1:12, header = FALSE)))
   # R2 <- as.vector(t(read.xlsx(workbook, sheetIndex = idx, rowIndex = 77:84, colIndex = 1:12, header = FALSE)))
   # R3 <- as.vector(t(read.xlsx(workbook, sheetIndex = idx, rowIndex = 86:93, colIndex = 1:12, header = FALSE)))

  
  # All substrates, raw values
  
   R1 <- as.numeric(as.vector(t(read.xlsx(workbook, sheetIndex = idx, rowIndex = 5:34, colIndex = 16:20, header = FALSE))))
   R2 <- as.numeric(as.vector(t(read.xlsx(workbook, sheetIndex = idx, rowIndex = 5:34, colIndex = 21:25, header = FALSE))))
   R3 <- as.numeric(as.vector(t(read.xlsx(workbook, sheetIndex = idx, rowIndex = 5:34, colIndex = 26:30, header = FALSE))))

  # Carbohydrates
  
  # R1 <- as.numeric(as.vector(t(read.xlsx(workbook, sheetIndex = idx, rowIndex = 5:34, colIndex = 46, header = FALSE))))
  # R2 <- as.numeric(as.vector(t(read.xlsx(workbook, sheetIndex = idx, rowIndex = 5:34, colIndex = 51, header = FALSE))))
  # R3 <- as.numeric(as.vector(t(read.xlsx(workbook, sheetIndex = idx, rowIndex = 5:34, colIndex = 56, header = FALSE))))
  
  # Polymers
   
   # R1 <- as.numeric(as.vector(t(read.xlsx(workbook, sheetIndex = idx, rowIndex = 5:34, colIndex = 47, header = FALSE))))
   # R2 <- as.numeric(as.vector(t(read.xlsx(workbook, sheetIndex = idx, rowIndex = 5:34, colIndex = 52, header = FALSE))))
   # R3 <- as.numeric(as.vector(t(read.xlsx(workbook, sheetIndex = idx, rowIndex = 5:34, colIndex = 57, header = FALSE))))
   
   # CARBOX & ACETIC
   
   # R1 <- as.numeric(as.vector(t(read.xlsx(workbook, sheetIndex = idx, rowIndex = 5:34, colIndex = 48, header = FALSE))))
   # R2 <- as.numeric(as.vector(t(read.xlsx(workbook, sheetIndex = idx, rowIndex = 5:34, colIndex = 53, header = FALSE))))
   # R3 <- as.numeric(as.vector(t(read.xlsx(workbook, sheetIndex = idx, rowIndex = 5:34, colIndex = 58, header = FALSE))))

   
   R1 <- R1[!is.na(R1)]
   R2 <- R2[!is.na(R2)]
   R3 <- R3[!is.na(R3)]
   

myDataSummary <- cbind(myDataSummary, R1)
colnames(myDataSummary)[length(colnames(myDataSummary))] <- paste0("W",timepoint,"-T",idx,"-R1")
myDataSummary <- cbind(myDataSummary, R2)
colnames(myDataSummary)[length(colnames(myDataSummary))] <- paste0("W",timepoint,"-T",idx,"-R2")
myDataSummary <- cbind(myDataSummary, R3)
colnames(myDataSummary)[length(colnames(myDataSummary))] <- paste0("W",timepoint,"-T",idx,"-R3")

# colnames(myDataSummary)[idx*3-2] <- paste0(idx,"-R1")
# colnames(myDataSummary)[idx*3-1] <- paste0(idx,"-R2")
# colnames(myDataSummary)[idx*3] <- paste0(idx,"-R3")

}

treatmentCategories <- NULL
for (idx in 1:number_of_sheets) {
  treatmentCategories <- c(treatmentCategories,idx,idx,idx)
  }

}


row.names(myDataSummary) <- compoundNames[,1]
myDataSummaryMetadata <- cbind.data.frame(myDataSummary, Substrate=gsub('[[:digit:]]+', '', compoundNames[,1]))


transformedDataSummary <- as.data.frame(t(myDataSummary))
transformedDataSummaryMeta <- as.data.frame(t(myDataSummaryMetadata))

CLPP_PCoA <- prcomp(scale(transformedDataSummary[,-1]))

###### FACTOEXTRA

fviz_pca_ind(CLPP_PCoA, col.ind = "cos2", # Color by the quality of representation
                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                        repel = TRUE     # Avoid text overlapping
             )

fviz_pca_ind(CLPP_PCoA, col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             label = "none"
)

fviz_pca_biplot(CLPP_PCoA, repel = TRUE,
                        col.var = "#2E9FDF", # Variables color
                        col.ind = "#696969"  # Individuals color
                )



####### GGPLOT2 #######


pca_plot_data1.2 <- as.data.frame(CLPP_PCoA$x[,1:2])
pca_plot_data2.3 <- as.data.frame(CLPP_PCoA$x[,2:3])

pca_plot_data1.2$Treatment <- as.factor(treatmentCategories)
pca_plot_data1.2$Week<-sapply(strsplit(row.names(pca_plot_data1.2), "-"), `[`, 1)
pca_plot_data2.3$Treatment <- as.factor(treatmentCategories)
pca_plot_data2.3$Week<-sapply(strsplit(row.names(pca_plot_data2.3), "-"), `[`, 1)

ggplot(pca_plot_data2.3, aes(x=PC2, y=PC3)) + geom_point(aes(color=Treatment), size=5)
ggplot(pca_plot_data1.2, aes(x=PC1, y=PC2)) + geom_point(aes(color=Treatment), size=5)
ggplot(pca_plot_data1.2, aes(x=PC1, y=PC2)) + geom_point(aes(color=Week), size=5)

library(phyloseq)

OTU <- otu_table(transformedDataSummary, taxa_are_rows = FALSE)

sampledata <- data.frame(Replicate= rep(1:3, times=7, each=1), Treatment=rep(c("1","2","3","4","5","6","7"), times=length(excelFiles), each=3), row.names=colnames(myDataSummary), stringsAsFactors=TRUE)
sampledata$Week <- sapply(strsplit(row.names(sampledata), "-"), `[`, 1)
sampledata$SampleID <- row.names(sampledata)
sampledata = sample_data(sampledata)
compoundTable<-data.frame(row.names=compoundNames[,1], Substrate=compoundNames[,1], Substrate_class=gsub('[[:digit:]]+', '', compoundNames[,1]))
compoundTable<-as.matrix(compoundTable)

physeq = phyloseq(otu_table(OTU), sample_data(sampledata), tax_table(compoundTable))
plot_bar(physeq, fill = "Substrate_class")

orderedByTreatment=colnames(myDataSummary)[c((1:3),(21+1:3),(42+1:3),(63+1:3),  3+(1:3),3+(21+1:3),3+(42+1:3),3+(63+1:3), 6+(1:3),6+(21+1:3),6+(42+1:3),6+(63+1:3), 9+(1:3),9+(21+1:3),9+(42+1:3),9+(63+1:3), 12+(1:3),12+(21+1:3),12+(42+1:3),12+(63+1:3), 15+(1:3),15+(21+1:3),15+(42+1:3),15+(63+1:3), 18+(1:3),18+(21+1:3),18+(42+1:3),18+(63+1:3))]

plot_heatmap(physeq, "NMDS", "bray", sample.order = orderedByTreatment, taxa.order= rev(row.names(myDataSummary)) )
plot_heatmap(physeq, "NMDS", "bray", sample.order = colnames(myDataSummary), taxa.order= rev(row.names(myDataSummary)) )

physeq.bray.dist <- distance(physeq, method="bray")
heatmap(data.matrix(physeq.bray.dist))

physeq.bray.dist <- distance(physeq, method="bray")
# heatmap(data.matrix(physeq.bray.dist))
physeq.ord <- ordinate(physeq, method="NMDS", distance="bray", formula=physeq.bray.dist)
adonisResults <- adonis2(physeq.bray.dist ~ Treatment, as(sample_data(physeq), "data.frame"))

plot_ordination(physeq, physeq.ord, color="Treatment") + geom_point(size=5)


plot_ordination(physeq, ordinate(physeq, "MDS", "bray"), color="Treatment") + geom_point(size=5)
plot_ordination(physeq, ordinate(physeq, "DCA", "bray"), color="Treatment") + geom_point(size=5)
plot_ordination(physeq, ordinate(physeq, "CCA", "bray"), color="Treatment") + geom_point(size=5)

plot_ordination(physeq, ordinate(physeq, "MDS", "bray"), color="Week", shape="Treatment") + geom_point(size=5)
plot_ordination(physeq, ordinate(physeq, "MDS", "bray"), color="Week") + geom_point(size=5)
plot_ordination(physeq, ordinate(physeq, "MDS", "bray"), color="Treatment", shape="Week") + geom_point(size=5)
 

## CCA



# physeq.ord.cca <- ordinate(physeq,"CCA", formula=physeq~Week)
# plot_ordination(physeq, physeq.ord.cca, type="split", color="Substrate_class") 
# ps_scores <- vegan::scores(physeq.ord.cca)
# sites <- data.frame(ps_scores$sites)
# sites$SampleID <- rownames(sites)
# sites <- sites %>% left_join(sample_data(physeq))
# 
# species <- data.frame(ps_scores$species)
# species$otu_id <- seq_along(colnames(otu_table(physeq)))
# tax <- tax_table(physeq)@.Data %>% data.frame(stringsAsFactors=FALSE)
# tax$otu_id <- seq_len(ncol(otu_table(physeq)))
# species <- species %>% left_join(tax)




cca.1 <- cca(t(myDataSummary) ~ myDataSummaryMetadata$Carb+myDataSummaryMetadata$Carb+myDataSummaryMetadata$Carb+myDataSummaryMetadata$Carb)
cca1.plot <- plot(cca.1,choices=c(1,2))
  

cca(bryceveg~elev+slope+av)





