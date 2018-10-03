# 
# Rgt1_10aaD.R
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
library(stringr)

samps <- c("3C_A1","3C_A2","3C_A3","3C_B1","3C_B2","3C_B3",
           "genomic_A1","genomic_A2","genomic_A3","genomic_B1","genomic_B2","genomic_B3")

# load data ----
combined <- read.table("Rgt1_10aaD_annotations.txt",col.names = "aa")
for (samp in samps) {
  temp <- read.table(paste(dir,"/Rgt1_10aaD_",samp,".counts.txt",sep=""),stringsAsFactors = FALSE, col.names=c("aa","count"))
  temp$norm <- temp$count/sum(temp$count)
  temp$count <- NULL
  colnames(temp) <- c("aa",samp)
  combined <- merge(combined, temp, by = "aa", all.x=TRUE)
}
combined[is.na(combined)] <- 0

# calculate ratios
combined$ratioA1 <- combined$`3C_A1`/combined$`genomic_A1`
combined$ratioA2 <- combined$`3C_A2`/combined$`genomic_A2`
combined$ratioA3 <- combined$`3C_A3`/combined$`genomic_A3`
combined$ratioB1 <- combined$`3C_B1`/combined$`genomic_B1`
combined$ratioB2 <- combined$`3C_B2`/combined$`genomic_B2`
combined$ratioB3 <- combined$`3C_B3`/combined$`genomic_B3`
combined$ratioA <- (combined$ratioA1 + combined$ratioA2 + combined$ratioA3)/3
combined$ratioB <- (combined$ratioB1 + combined$ratioB2 + combined$ratioB3)/3
combined$ratio <- (combined$ratioA + combined$ratioB)/2

# load predicted disorder data
idr <- read.table("Rgt1_IUPred2.txt",comment.char = "#",col.names=c("pos","aa","IUPred","Anchor"))
idr$aa10 <- floor((idr$pos + 9)/10)*10

# load secondary structure predictions
jpred0 <- read.table("Rgt1_1-400_Jpred4/jp__4c8Orq.txt",sep = "\t",header=F, col.names=c("jnetpred"))
jpred0$aa <- as.numeric(rownames(jpred0))
jpred0$jnetpred <- factor(jpred0$jnetpred,levels=c("-","H","E"))

jpred <- read.table("Rgt1_401-1170_Jpred4/jp_PsBbSbJ.txt",sep = "\t",header=F, col.names=c("jnetpred"))
jpred$aa <- as.numeric(rownames(jpred)) + 400
jpred$jnetpred <- factor(jpred$jnetpred,levels=c("-","H","E"))

jpredall <- rbind(jpred0,jpred)

# generate figures ----
brewercols <- brewer.pal(6,"Set1")

# for paper, Figure 4D
pdf(paste(out,"/Rgt1_10aaD_regulatory_domain_overlay.pdf",sep=""),3.8,1.8)
ggplot(combined) + 
  geom_rect(aes(xmin=520,xmax=550,ymin=-Inf,ymax=Inf), fill="tan1") +
  geom_rect(aes(xmin=640,xmax=720,ymin=-Inf,ymax=Inf), fill="tan1") +
  geom_rect(aes(xmin=750,xmax=800,ymin=-Inf,ymax=Inf), fill="tan1") +
  geom_rect(aes(xmin=810,xmax=830,ymin=-Inf,ymax=Inf), fill="tan1") +
  geom_bar(aes(x=aa+5,y=log2(ratio)),stat="identity",fill=brewercols[1],width=10) + theme_classic() + 
  geom_segment(aes(x=aa+5,xend=aa+5,y=log2(ratioA),yend=log2(ratioB)),size=0.2) + 
  geom_hline(aes(yintercept=0)) +
  xlab("Rgt1 amino acid position") + ylab("Log2 3C/Genomic") + scale_x_continuous(expand=c(0,0),limits=c(0,1170),breaks=seq(0,1170,250)) + paper.font
dev.off()

# for paper, Figure 4--figure supplement 1
pdf(paste(out,"/Rgt1_10aaD_iupred_jnetpred_overlay.pdf",sep=""),6,3)
ggplot(combined) + 
  geom_rect(data=jpredall,aes(xmin=aa-1,xmax=aa,ymin=-Inf,ymax=Inf,fill=jnetpred)) + scale_fill_manual(values=c("white","pink","lightblue")) + 
  geom_bar(aes(x=aa+5,y=log2(ratio)),stat="identity",width=10,fill=brewercols[1]) + 
  geom_segment(aes(x=aa+5,xend=aa+5,y=log2(ratioA),yend=log2(ratioB)),size=0.2) + 
  theme_classic() + xlab("Rgt1 amino acid position") + ylab("Log2 (3C/Genomic)") + scale_x_continuous(expand=c(0,0),limits=c(0,1170),breaks=seq(0,1170,250)) + 
  geom_line(data=idr,aes(x=pos,y=IUPred*4-3),color="black") + 
  theme(legend.position="none",text=element_text(size=10),axis.text=element_text(size=10,color="black")) + 
  scale_y_continuous(sec.axis = sec_axis(~./4+.75,name="Predicted disorder"))
dev.off()

# no overlay
pdf(paste(out,"/Rgt1_10aaD.pdf",sep=""),3.8,1.8)
ggplot(combined) + 
  geom_bar(aes(x=aa+5,y=log2(ratio)),stat="identity",fill=brewercols[1],width=10) + theme_classic() + 
  geom_segment(aes(x=aa+5,xend=aa+5,y=log2(ratioA),yend=log2(ratioB)),size=0.2) + 
  geom_hline(aes(yintercept=0)) +
  xlab("Rgt1 amino acid position") + ylab("Log2 3C/Genomic") + scale_x_continuous(expand=c(0,0),limits=c(0,1170),breaks=seq(0,1170,250)) + paper.font
dev.off()

# all annotations from Polish et al, Genetics 2005
pdf(paste(out,"/Rgt1_10aaD_Polish2005_overlay.pdf",sep=""),3.8,1.8)
ggplot(combined) + 
  geom_rect(aes(xmin=150,xmax=160,ymin=-Inf,ymax=Inf), fill="pink") +
  geom_rect(aes(xmin=210,xmax=250,ymin=-Inf,ymax=Inf), fill="pink") +
  geom_rect(aes(xmin=310,xmax=320,ymin=-Inf,ymax=Inf), fill="lightgreen") +
  geom_rect(aes(xmin=320,xmax=360,ymin=-Inf,ymax=Inf), fill="lightblue") +
  geom_rect(aes(xmin=370,xmax=380,ymin=-Inf,ymax=Inf), fill="lightblue") +
  geom_rect(aes(xmin=400,xmax=410,ymin=-Inf,ymax=Inf), fill="lightgreen") +
  geom_rect(aes(xmin=520,xmax=550,ymin=-Inf,ymax=Inf), fill="tan1") +
  geom_rect(aes(xmin=640,xmax=720,ymin=-Inf,ymax=Inf), fill="tan1") +
  geom_rect(aes(xmin=750,xmax=800,ymin=-Inf,ymax=Inf), fill="tan1") +
  geom_rect(aes(xmin=810,xmax=830,ymin=-Inf,ymax=Inf), fill="tan1") +
  geom_rect(aes(xmin=330,xmax=340,ymin=-Inf,ymax=Inf), fill="lightgrey") +
  geom_rect(aes(xmin=440,xmax=460,ymin=-Inf,ymax=Inf), fill="lightgrey") +
  geom_rect(aes(xmin=730,xmax=750,ymin=-Inf,ymax=Inf), fill="lightgrey") +
  geom_rect(aes(xmin=800,xmax=810,ymin=-Inf,ymax=Inf), fill="lightgrey") +
  geom_rect(aes(xmin=880,xmax=900,ymin=-Inf,ymax=Inf), fill="lightgrey") +
  geom_rect(aes(xmin=920,xmax=970,ymin=-Inf,ymax=Inf), fill="lightgrey") +
  geom_rect(aes(xmin=1000,xmax=1040,ymin=-Inf,ymax=Inf), fill="lightgrey") +
  geom_bar(aes(x=aa+5,y=log2(ratio)),stat="identity",fill=brewercols[1],width=10) + theme_classic() + 
  geom_segment(aes(x=aa+5,xend=aa+5,y=log2(ratioA),yend=log2(ratioB)),size=0.2) + 
  geom_hline(aes(yintercept=0)) +
  xlab("Rgt1 amino acid position") + ylab("Log2 3C/Genomic") + scale_x_continuous(expand=c(0,0),limits=c(0,1170),breaks=seq(0,1170,250)) + paper.font
dev.off()
