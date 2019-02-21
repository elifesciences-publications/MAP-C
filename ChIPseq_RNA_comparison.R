# 
# ChIPseq_RNA_comparison.R
# 
# Seungsoo Kim
# February 20, 2019

# directories ----
setwd("/Volumes/shendure-vol8/projects/mutagenesis.3C/nobackup/MAP-C")
dir <- "nobackup/ChIPseq/macs"
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

fig <- theme_classic() + theme(text=element_text(size=8),
                               axis.text=element_text(size=8,color="black"),
                               legend.text=element_text(size=8), strip.text=element_text(size=8), strip.background = element_blank())

tfs <- c("Leu3","Sdd4","Rgt1")
TFs <- c("LEU3","YPR022C","RGT1")

for (tfi in 1:length(tfs)) {
  tf <- tfs[tfi]
  TF <- TFs[tfi]
  
  # load saturated peaks labeled by downstream gene
  x.sat <- read.table(paste(dir,"/",tf,"_saturated_merged_peaks_labeled.bed",sep=""),col.names=c("chr","start","end","fe","Ldir","Lsgd","Lname","Ldist","Rdir","Rsgd","Rname","Rdist","name"))
  x.sat <- x.sat[order(-x.sat$fe),]
  
  # load exponential peaks labeled by downstream gene
  x.exp <- read.table(paste(dir,"/",tf,"_exponential_merged_peaks_labeled.bed",sep=""),col.names=c("chr","start","end","fe","Ldir","Lsgd","Lname","Ldist","Rdir","Rsgd","Rname","Rdist","name"))
  x.exp <- x.exp[order(-x.exp$fe),]
  
  # load RNA-seq data
  y <- read.table(paste("nobackup/",tf,"D_RNAfc.txt",sep=""))
  
  # merge saturated datasets
  m.sat <- merge(y,x.sat,by.x="name",by.y="name",all.x=T)
  m.sat[is.na(m.sat$fe),]$fe <- 0.5
  
  # merge exponential datasets
  m.exp <- merge(y,x.exp,by.x="name",by.y="name",all.x=T)
  m.exp[is.na(m.exp$fe),]$fe <- 0.5
  
  # for paper, Figure 5--figure supplement 2A
  pdf(paste("nobackup/2019-01-20_ChIPseq/figures/",tf,"_saturated_RNA_comparison_fc.pdf",sep=""),2.2,2.2)
  print(ggplot(subset(m.sat,name != TF)) + geom_point(aes(x=log2(fe),y=log2FoldChange),size=0.3) + fig + 
          xlab(expression(paste("ChIP ",Log[2]," Fold Enrichment"))) + ylab(expression(paste("RNA ",Log[2]," Fold Change"))) + theme(legend.position="none"))
  dev.off()
  
  # for paper, Figure 5--figure supplement 2B
  pdf(paste("nobackup/2019-01-20_ChIPseq/figures/",tf,"_exponential_RNA_comparison_fc.pdf",sep=""),2.2,2.2)
  print(ggplot(subset(m.exp,name != TF)) + geom_point(aes(x=log2(fe),y=VALUE),size=0.3) + fig + 
          xlab(expression(paste("ChIP ",Log[2]," Fold Enrichment"))) + ylab(expression(paste("RNA ",Log[2]," Fold Change"))) + theme(legend.position="none"))
  dev.off()
    
}



