# 
# fixed-locus_ko.R
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

samps <- c("3C_1","3C_2","3C_3","offtarget_1","offtarget_2","offtarget_3","genomic_1","genomic_2","genomic_3")
classes <- c("LEU3","RGT1","MOT3","TF","Nup","Neutral","WT","blank")
paper.font <- theme(text=element_text(size=8), axis.text=element_text(size=8,color="black"), legend.text=element_text(size=8))

# load data ----
combined <- read.table("fixed-locus_ko_annotations.txt", col.names=c("name","systname","gene","class","growth"),stringsAsFactors = FALSE)
combined$class <- factor(combined$class,levels=classes)
combined <- subset(combined,growth=="y")
for (samp in samps) {
  temp <- read.table(paste(dir,"/fixed-locus_ko_",samp,".counts.txt",sep=""),stringsAsFactors = FALSE,col.names = c("name","systname","count"))
  temp$norm <- temp$count/sum(temp$count)
  temp$count <- NULL
  temp$systname <- NULL
  temp <- subset(temp,norm > 0)
  colnames(temp) <- c("name",samp)
  combined <- merge(combined, temp, by = "name", all.x=TRUE)
}
#combined[is.na(combined)] <- 0

# calculate ratios
combined$ratio1 <- combined$`3C_1`/combined$`genomic_1`
combined$ratio2 <- combined$`3C_2`/combined$`genomic_2`
combined$ratio3 <- combined$`3C_3`/combined$`genomic_3`
combined$ratio <- (combined$ratio1 + combined$ratio2 + combined$ratio3)/3
combined$off1 <- combined$`offtarget_1`/combined$`genomic_1`
combined$off2 <- combined$`offtarget_2`/combined$`genomic_2`
combined$off3 <- combined$`offtarget_3`/combined$`genomic_3`
combined$ratiosd <- apply(as.matrix(cbind(combined$ratio1,combined$ratio2,combined$ratio3)),1,sd)
combined$inputsd <- apply(as.matrix(cbind(combined$`genomic_1`,combined$`genomic_2`,combined$`genomic_3`)),1,sd)

# generate figures ----
brewercols <- brewer.pal(6,"Set1")
classcols <- c(brewercols[2],brewercols[5],brewercols[1],"lightgrey",brewercols[4],"grey40","black","white")

# for paper, Figure 2B
pdf(paste(out,"/fixed-locus_ko_pairratio_hist_filtered.pdf",sep=""),3.2,2)
ggplot(subset(combined,(`genomic_1`+`genomic_2`+`genomic_3`)/3 > .003)) + 
  geom_histogram(aes(x=log2(ratio),fill=class),binwidth=.1) + 
  scale_fill_manual(values=classcols,name="") + theme_classic() + 
  scale_y_continuous(expand=c(0,0)) + xlab(expression(paste(Log[2]," 3C/Genomic"))) + ylab("Barcoded strains") + paper.font +
  theme(legend.key.size = unit(0.3,"cm"), legend.position = c(0.15,0.8))
dev.off()

# for paper, Figure 2--figure supplement 2B
# input v pairing ratio scatter
#growing <- subset(combined,growth=="y")
xshift <- -0.0002
highlight <- c("VHR1","CBF1")
pdf(paste(out,"/fixed-locus_ko_input_v_pairratio_scatter_errorsd.pdf",sep=""),7,3)
ggplot(subset(combined,growth=="y")[order(subset(combined,growth=="y")$class),]) + 
  geom_segment(aes(x=(`genomic_1`+`genomic_2`+`genomic_3`)/3-inputsd,xend=(`genomic_1`+`genomic_2`+`genomic_3`)/3+inputsd,
                   y=log2(ratio),yend=log2(ratio),color=class)) +
  geom_segment(aes(x=(`genomic_1`+`genomic_2`+`genomic_3`)/3,xend=(`genomic_1`+`genomic_2`+`genomic_3`)/3,
                   y=log2(ratio-ratiosd),yend=log2(ratio+ratiosd),color=class)) +
  geom_point(aes(x=(`genomic_1`+`genomic_2`+`genomic_3`)/3,y=log2(ratio),color=class)) + 
  scale_color_manual(values=classcols, name="") + theme_bw() + xlab("Average genomic frequency") + ylab(expression(paste(Log[2]," 3C/Genomic"))) + paper.font +
  geom_text(data=subset(combined,gene %in% highlight),aes(x=(`genomic_1`+`genomic_2`+`genomic_3`)/3+xshift,y=log2(ratio),label=gene),hjust=1,size=8*5/14)
dev.off()

