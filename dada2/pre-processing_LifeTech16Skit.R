# Plots the aligned locations of 16S amplicons from the LifeTech kit
args<-commandArgs(TRUE)
# library(ggplot2)

alignment <- args[1]

data <- read.table(alignment, sep="\t", header=T)
data$Number<-c(1:nrow(data))

data$Region <- ifelse(data$TemplateEnd<=370,"V2","")
data$Region <- ifelse(data$TemplateStart>200 & data$TemplateEnd<600,"V3",data$Region)
data$Region <- ifelse(data$TemplateStart>450 & data$TemplateEnd<850,"V4",data$Region)
data$Region <- ifelse(data$TemplateStart>850 & data$TemplateEnd<1200,"V6-7",data$Region)
data$Region <- ifelse(data$TemplateStart>1000 & data$TemplateEnd<1370,"V8",data$Region)
data$Region <- ifelse(data$TemplateStart>1200,"V9",data$Region)
data$Region <- ifelse(data$TemplateStart>1450,NA,data$Region)

newFilename <- gsub(".align.report","",alignment)

write.table(data, file=paste0(newFilename,".regions"), sep="\t", quote=FALSE)

# dataTop <- head(data, n=500)

# ggplot(dataTop) + geom_segment(aes(x=TemplateStart,xend=TemplateEnd,y=Number,yend=Number, color=Region),alpha=0.5) + scale_x_continuous(breaks=seq(0, 1500, 100), limits=c(0,1500))
# ggplot(data) + geom_segment(aes(x=TemplateStart,xend=TemplateEnd,y=Number,yend=Number, color=Region),alpha=0.5) + scale_x_continuous(breaks=seq(0, 1500, 100), limits=c(0,1500))

