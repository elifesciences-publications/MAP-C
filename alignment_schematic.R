# 
# alignment_schematic.R
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
library("GGally")
library(stringr)

brewercols <- brewer.pal(6,"Set1")
classcols <- c(brewercols[2],brewercols[3],brewercols[5],"gray50")

specs <- c("Scer","Spar","Smik","Skud","Suva")
TFs <- c("Leu3","Sdd4","Rgt1")

x <- read.table("alignment.txt",fill=T,header=F)
colnames(x) <- c("specno","start","end","tf","score")
x$species <- x$specno
x$specno <- as.numeric(factor(x$specno,levels=specs))

species <- c("cerevisiae","paradoxus","mikatae","kudriavzevii","uvarum")
pdf(paste(out,"/alignment_schematic.pdf",sep=""),6.5,2)
ggplot() + geom_segment(data=subset(x,tf=="full"),aes(x=start,xend=end,y=specno,yend=specno)) + 
  geom_rect(data=subset(x,tf!="full"),aes(xmin=start,xmax=end,ymin=specno,ymax=specno-0.5*score/100,fill=factor(tf,levels=TFs),alpha=score/100),color="black") + 
  scale_fill_manual(values=classcols) + scale_y_reverse(breaks=1:5,labels=species) + labs(alpha="Score",fill="TF") +
  theme_classic() + theme(text=element_text(size=8),axis.line.y = element_blank(),axis.ticks.y = element_blank(), axis.text.y=element_text(size=8,face="italic",color="black"), axis.text.x=element_text(size=8,color="black"),legend.key.size = unit(0.3,"cm")) + ylab("") + xlab("Position (bp)")
dev.off()
