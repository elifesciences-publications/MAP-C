# 
# Rgt1_muts_schematic.R
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

paper.font <- theme(text=element_text(size=8), axis.text=element_text(size=8,color="black"), legend.text=element_text(size=8))

domains <- read.table("Rgt1_domains.txt", sep="\t", col.names = c("domain","start","end"))
muts <- read.table("Rgt1_mutations.txt", sep="\t", col.names = c("mut","start","end","offset"))
muts$mut <- factor(muts$mut,levels=c("Zn finger", "Q/N-rich","C-terminal","S88A","S758A"))

brewercols <- brewer.pal(6,"Set1")
oranges <- brewer.pal(7,"Oranges")
cols <- oranges[c(2,3,5,6,7)]

# for Figure 4E
pdf(paste(out,"/Rgt1_muts_schematic.pdf",sep=""),3.2,1)
ggplot() + 
  geom_segment(data=subset(domains,domain=="full"),aes(x=start,xend=end,y=0,yend=0)) + 
  geom_rect(data=muts,aes(xmin=start-0.5,xmax=end+0.5,ymin=-0.03,ymax=0.03+(offset-1)*0.025,fill=mut),color="black") + 
  geom_text(data=muts,aes(label=mut,x=(start+end)/2,y=0.05*offset+0.055), size = 8*5/14) + theme_classic() + 
  scale_fill_manual(values=cols,name="") + paper.font +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position="none") + 
  ylab("") + xlab("Rgt1 amino acid position") + ylim(-0.05,0.2) +
  theme(plot.margin=unit(c(0,0.3,0,0),"cm"))
dev.off()