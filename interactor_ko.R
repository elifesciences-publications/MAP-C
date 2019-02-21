# 
# interactors_ko.R
# 
# Seungsoo Kim
# February 19, 2019

# directories ----
setwd("/Volumes/shendure-vol8/projects/mutagenesis.3C/nobackup/MAP-C")
dir <- "data"
out <- "figures"

# load libraries ----
library(ggplot2)
library(dplyr)
library(reshape)
library(RColorBrewer)
library(data.table)
library(grid)
library(gplots)
library(scales)
library(permute)
library(gridExtra)
library(stringr)

samps <- c("3C_1","3C_2","3C_3","genomic_1","genomic_2","genomic_3")
classes <- c("RNQ1","TUP1","SSN6","MTH1","STD1","LEU3","RGT1","SDD4","MOT3","Neutral","WT")
paper.font <- theme(text=element_text(size=8), axis.text=element_text(size=8,color="black"), legend.text=element_text(size=8))

# load data ----
combined <- read.table("interactor_ko_annotations.txt", col.names=c("name","systname","gene","class","growth"))
combined$class <- factor(combined$class,levels=classes)
for (samp in samps) {
  temp <- read.table(paste(dir,"/interactor_ko_",samp,".counts.txt",sep=""),stringsAsFactors = FALSE,col.names = c("name","systname","count"))
  temp$norm <- temp$count/sum(temp$count)
  temp$count <- NULL
  temp$systname <- NULL
  colnames(temp) <- c("name",samp)
  combined <- merge(combined, temp, by = "name", all=TRUE)
}
combined[is.na(combined)] <- 0

# calculate ratios
combined$ratio1 <- combined$`3C_1`/combined$`genomic_1`
combined$ratio2 <- combined$`3C_2`/combined$`genomic_2`
combined$ratio3 <- combined$`3C_3`/combined$`genomic_3`
combined$ratio <- (combined$ratio1 + combined$ratio2 + combined$ratio3)/3
combined$ratiosd <- apply(as.matrix(cbind(combined$ratio1,combined$ratio2,combined$ratio3)),1,sd)
combined$inputsd <- apply(as.matrix(cbind(combined$`genomic_1`,combined$`genomic_2`,combined$`genomic_3`)),1,sd)

# generate figures ----
brewercols <- brewer.pal(6,"Set1")
purples <- brewer.pal(5,"Purples")
classcols <- c(purples,brewercols[2],brewercols[5],brewercols[3],brewercols[1],"grey40","black")

# for paper, Figure 4C
pdf(paste(out,"/interactor_ko_bar_scatter.pdf",sep=""),3.2,2.1)
ggplot(subset(combined,(`genomic_1` + `genomic_2` + `genomic_3` > 0.001) & (class %in% classes))) + 
  geom_bar(aes(x=class,y=log2(ratio),fill=class),stat="summary",fun.y="median",color="black") + 
  geom_point(aes(x=class,y=log2(ratio)),position="jitter",size=0.5) + 
  scale_fill_manual(values=classcols) + theme_classic() + theme(legend.position="none") + xlab("") + ylab(expression(paste(Log[2]," 3C/Genomic"))) + paper.font +
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=0), axis.text.y=element_text(hjust=0))
dev.off()
