# 
# in-gene_ko.R
# 
# Seungsoo Kim
# September 26, 2018

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

samps <- c("genomic_1","genomic_2","offtarget_1","offtarget_2","3C_1","3C_2")
brewercols <- brewer.pal(n=6,name="Set1")
paper.font <- theme(text=element_text(size=8), axis.text=element_text(size=8,color="black"), legend.text=element_text(size=8))

cendist <- read.table("in-gene_ko_cendist.txt", col.names=c("gene","cendist"))
groups <- read.table("in-gene_ko_groups.txt",stringsAsFactors = FALSE, col.names=c("gene","group"))

# load data ----
combined <- merge(groups,cendist,by="gene",all=TRUE)
for (samp in samps) {
  temp <- read.table(paste(dir,"/in-gene_ko_",samp,".counts.txt",sep=""),stringsAsFactors = FALSE, col.names=c("gene","count"))
  temp$norm <- temp$count/sum(temp$count)
  temp$count <- NULL
  colnames(temp) <- c("gene",samp)
  combined <- merge(combined, temp, by = "gene", all=TRUE)
}
combined[is.na(combined)] <- 0

# calculate ratios
combined$ratio1 <- combined$`3C_1`/combined$`genomic_1`
combined$ratio2 <- combined$`3C_2`/combined$`genomic_2`

# generate figures ----
combined$groupfac <- factor(combined$group)
levels(combined$groupfac) <- c("YDR207C(UME6)","YER028C(MIG3)","YGL035C(MIG1)","YGL209W(MIG2)","YKL038W(RGT1)","YKL062W(MSN4)","YLR451W(LEU3)","YMR037C(MSN2)","YOR113W(AZF1)","YPL248C(GAL4)")

# for paper, Figure 2--figure supplement 1B
pdf(paste(out,"/in-gene_ko_effect_of_centromere_distance.pdf",sep=""),2.5,4)
ggplot(subset(combined,`genomic_1` > .002 & `genomic_2` > .002)) + 
  geom_point(aes(x=cendist,y=ratio1,color=gene!=group)) + 
  geom_point(aes(x=cendist,y=ratio2,color=gene!=group)) + 
  geom_segment(aes(x=cendist,xend=cendist,y=ratio1,yend=ratio2,color=gene!=group)) + 
  theme_classic() + scale_color_manual(values=c("red","black"),labels=c("TF","Neighboring gene"),name="") + 
  xlab("Distance from centromere (bp)") + ylab("3C/Genomic") + 
  theme(text=element_text(size=8,color="black"),legend.position=c(0.8,0.2))
dev.off()

# for paper, Figure 2--figure supplement 1C
pdf(paste(out,"/in-gene_ko_pairing_ratio_v_input_bygroup.pdf",sep=""),6.6,3.5)
ggplot(combined[order(combined$gene!=combined$group),]) + 
  geom_point(aes(x=`genomic_1`,y=ratio1,color=gene!=group)) + 
  geom_point(aes(x=`genomic_2`,y=ratio2,color=gene!=group)) + 
  geom_segment(aes(x=`genomic_1`,xend=`genomic_2`,y=ratio1,yend=ratio2,color=gene!=group)) + 
  theme_bw() + scale_x_log10(limits = c(-.004,.071)) + 
  facet_wrap(~groupfac,ncol=5) + 
  scale_color_manual(values=c("red","black"),labels=c("TF","Neighboring gene"),name="") + 
  ylab("3C/Genomic") + xlab("Frequency in genomic") + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),                                                                                                                                                                                                                                       labels = trans_format("log10", math_format(10^.x))) + annotation_logticks(sides="l") + theme(text=element_text(size=8,color="black"),strip.background=element_blank())
dev.off()

