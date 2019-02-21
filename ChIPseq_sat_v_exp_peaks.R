# 
# ChIPseq_sat_v_exp_peaks.R
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

tfs <- c("Leu3","Sdd4","Rgt1")
fig <- theme_classic() + theme(text=element_text(size=8),
                               axis.text=element_text(size=8,color="black"),
                               legend.text=element_text(size=8), strip.text=element_text(size=8), strip.background = element_blank())

for (tf in tfs) {
    x <- read.table(paste("nobackup/ChIPseq/macs/",tf,"_saturated_merged_overlap_exponential_peaks.bed",sep=""),col.names=c("chr","start","end","fe","olap"))
    y <- read.table(paste("nobackup/ChIPseq/macs/",tf,"_exponential_merged_overlap_saturated_peaks.bed",sep=""),col.names=c("chr","start","end","fe","olap"))
    y2 <- subset(y,olap==0)
    y2$olap <- y2$fe
    y2$fe <- 0
    z <- rbind(x,y2)

    # for Figure 2--figure supplement 3A
    png(paste("nobackup/2019-01-20_ChIPseq/figures/",tf,"_sat_v_exp_peaks_scatter.png",sep=""),2.1,2.1,units="in",res=1200)
    print(ggplot(z) + geom_point(data=subset(z,olap!=0 & fe!=0),aes(x=fe,y=olap,color=chr=="chrXIII" & start > 850000 & end < 855000),size=0.3) + 
            geom_point(data=subset(z,fe==0),aes(x=0.5,y=olap),size=0.3,position=position_jitter(height = 0,width=0.4)) +
            geom_point(data=subset(z,olap==0),aes(x=fe,y=0.5,color=chr=="chrXIII" & start > 850000 & end < 855000),size=0.3,position=position_jitter(height = 0.4,width=0)) + 
            geom_hline(aes(yintercept=1)) + geom_vline(aes(xintercept=1)) + scale_color_manual(values=c("black","red")) +
            fig + theme(legend.position="none") +
            xlab("Saturated Fold Enrichment") + ylab("Exponential Fold Enrichment"))
    dev.off()
}


