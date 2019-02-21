# 
# ChIPseq_motifs.R
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

for (tf in c("Leu3","Sdd4","Rgt1")) {
  for (cond in c("saturated","exponential")) {
    for (pval in c("0.001","0.0005")) {
      # load data - table of ChIP-seq peaks, including counts of motifs
      x <- read.table(paste(dir,"/",tf,"_",cond,"_merged_summits_motifs_",pval,".bed",sep=""),col.names=c("chr","start","end","fe","motifs"))
      
      # for paper, Figure 2--figure supplement 3B and C
      # boxplots with overlaid scatter plots
      png(paste(out,"/ChIPseq_",tf,"_",cond,"_merged_motifs_",pval,"_boxscatter.png",sep=""),2.1,2.1,units="in",res=1200)
      print(ggplot(x) + geom_boxplot(aes(x=motifs,y=log10(fe),group=motifs),outlier.shape = NA) + 
              geom_point(aes(x=motifs,y=log10(fe)),position=position_jitter(0.3),size=0.1) + fig + 
              xlab("Number of motifs") + ylab(expression(paste(Log[10]," Fold Enrichment"))))
      dev.off()
    }
    
  }
}
