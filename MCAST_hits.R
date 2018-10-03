# 
# MCAST_hits.R
# 
# Seungsoo Kim
# September 26, 2018

# directories ----
setwd("/Volumes/shendure-vol8/projects/mutagenesis.3C/nobackup/MAP-C")
dir <- "data"
out <- "figures"

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

brewercols <- brewer.pal(n=6,name="Set1")
classcols <- c(brewercols[2],brewercols[3],brewercols[5])
paper.font <- theme(text=element_text(size=8), axis.text=element_text(size=8,color="black"), legend.text=element_text(size=8))

TFs <- c("Leu3","Sdd4","Rgt1")
x <- read.table("Scer_0.0005_50_Leu3_Sdd4_Rgt1_LR.motifs.txt")
colnames(x) <- c("chr","chrstart","chrend","locus","score","motifno","motifstart","motifend","motifscore","motifdir","TF","pval")
x <- x[order(x$pval),]
x$TF <- factor(x$TF,levels=TFs)
x$genome <- "Scer"

y <- read.table("Sbay_0.0005_50_Leu3_Sdd4_Rgt1_LR.motifs.txt")
colnames(y) <- c("chr","chrstart","chrend","locus","score","motifno","motifstart","motifend","motifscore","motifdir","TF","pval")
y <- y[order(y$pval),]
y$TF <- factor(y$TF,levels=TFs)
y$genome <- "Suva"

z <- rbind(x,y)

pdf(paste(out,"/MCAST_hits.pdf",sep=""),2.5,1.7)
ggplot(z) + geom_rect(aes(xmin=motifstart, xmax=motifend+1, ymin=0, ymax=-log10(pval), fill=TF), color="black") + 
  geom_segment(aes(x=0,xend=chrend-chrstart+1,y=0,yend=0)) +
  facet_grid(motifno ~ genome) + theme_classic() + paper.font +
  theme(panel.spacing = unit(0.35,"cm"), strip.background = element_blank(), strip.text = element_blank(), 
        legend.position = c(0.75,0.35), legend.background = element_blank(), legend.key.size=unit(0.3,"cm")) + 
  xlab("") + ylab("") + 
  scale_fill_manual(values=classcols,name="") + scale_y_continuous(breaks=c(0,3,6))
dev.off()
