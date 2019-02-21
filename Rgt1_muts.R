# 
# Rgt1_muts.R
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

samps <- c("genomic_1","genomic_2","genomic_3","3C_1","3C_2","3C_3")

classes <- c("RGT1_WT","RGT1_ZF","RGT1_QN","RGT1_CT","RGT1_S88A","RGT1_S758A","RGT1","LEU3","SDD4","MOT3","Neutral","WT")

paper.font <- theme(text=element_text(size=8), axis.text=element_text(size=8,color="black"), legend.text=element_text(size=8))

# load data ----
combined <- read.table("Rgt1_muts_annotations.txt", col.names=c("name","systname","gene","class","growth"))
combined$class <- factor(combined$class,levels=classes)
for (samp in samps) {
  temp <- read.table(paste(dir,"/Rgt1_muts_",samp,".counts.txt",sep=""),stringsAsFactors = FALSE,col.names = c("name","systname","count"))
  temp$norm <- temp$count/sum(temp$count)
  temp$count <- NULL
  temp$systname <- NULL
  colnames(temp) <- c("name",samp)
  combined <- merge(combined, temp, by = "name", all=TRUE)
}

# calculate ratios
combined$ratio1 <- combined$`3C_1`/combined$`genomic_1`
combined$ratio2 <- combined$`3C_2`/combined$`genomic_2`
combined$ratio3 <- combined$`3C_3`/combined$`genomic_3`
combined$ratio <- (combined$ratio1 + combined$ratio2 + combined$ratio3)/3
combined$ratiosd <- apply(as.matrix(cbind(combined$ratio1,combined$ratio2,combined$ratio3)),1,sd)
combined$inputsd <- apply(as.matrix(cbind(combined$`genomic_1`,combined$`genomic_2`,combined$`genomic_3`)),1,sd)

# generate figures ----
brewercols <- brewer.pal(6,"Set1")
oranges <- brewer.pal(7,"Oranges")
classcols <- c(oranges[1:3],oranges[5:7],brewercols[5],brewercols[2],brewercols[3],brewercols[1],"grey40","black")

# for paper, Figure 4F
pdf(paste(out,"/Rgt1_muts_bar_scatter.pdf",sep=""),3.2,2.1)
ggplot(subset(combined,(`genomic_1` + `genomic_2` + `genomic_3` > 0.001) & (class %in% classes))) + 
  geom_bar(aes(x=class,y=log2(ratio),fill=class),stat="summary",fun.y="median",color="black") + 
  geom_point(aes(x=class,y=log2(ratio)),position="jitter",size=0.5) + 
  scale_fill_manual(values=classcols) + theme_classic() + theme(legend.position="none") + xlab("") + ylab(expression(paste(Log[2]," 3C/Genomic"))) + paper.font +
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=0), axis.text.y=element_text(hjust=0))
dev.off()
