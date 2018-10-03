# 
# subsequences.R
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

# load data ----
samps <- c("3C_1","3C_2","genomic_1","genomic_2","offtarget_1","offtarget_2")
libs <- c("3C","offtarget","genomic")
paper.font <- theme(text=element_text(size=8), axis.text=element_text(size=8,color="black"), legend.text=element_text(size=8))

combined.cov <- data.frame(pos=integer(), cov=integer(), samp=character())
for (samp in samps) {
  cov <- read.table(paste(dir,"/subsequences_",samp,".cov.txt",sep=""), col.names=c("seq","pos","cov"))
  cov <- cov[,-1]
  cov$cov <- cov$cov/sum(cov$cov) # normalize coverage by total
  cov$samp <- samp
  combined.cov <- rbind(combined.cov,cov)
}
combined.cov$lib <- factor(str_sub(combined.cov$samp,1,str_length(combined.cov$samp)-2),levels=libs)

# generate figures ----  
paircols <- brewer.pal(6,"Paired")
linecols <- c(paircols[6],"grey50")
linelabs <- c("3C","Genomic")

# for paper, Figure 1B
# coverage of 3C and genomic libraries
pdf(paste(out,"/subsequences_coverage.pdf",sep=""),3.3,2)
ggplot(subset(combined.cov,lib!="offtarget")) + 
  geom_line(aes(x=pos,y=cov,group=samp,color=lib)) + 
  theme_classic() + scale_color_manual(values=linecols,labels=linelabs,name="") + 
  xlab("Position (bp)") + ylab("Relative coverage") + paper.font + 
  theme(legend.key.height=unit(1,"line"),legend.position=c(0.85,0.8)) + ylim(0,0.0053)
dev.off()

# for paper, Figure 1--figure supplement 1B
# coverage of offtarget and genomic libraries
pdf(paste(out,"/subsequences_coverage_offtarget.pdf",sep=""),3.3,2)
ggplot(subset(combined.cov,lib!="3C")) + 
  geom_line(aes(x=pos,y=cov,group=samp,color=lib)) + 
  theme_classic() + scale_color_manual(values=c(paircols[2],"grey50"),labels=linelabs,name="") + 
  xlab("Position (bp)") + ylab("Relative coverage") + paper.font + 
  theme(legend.key.height=unit(1,"line"),legend.position=c(0.85,0.8)) + ylim(0,0.0053)
dev.off()

# for paper, Figure 1--figure supplement 1C
# coverage of 3C and genomic libraries overlaid with highlights of regions used for ectopic pairing and nonpairing sequences
pdf(paste(out,"/subsequences_coverage_ectopic.pdf",sep=""),3.3,2)
ggplot(subset(combined.cov,lib!="offtarget"))+ 
  geom_rect(aes(xmin=509.5,xmax=695.5,ymin=-Inf,ymax=Inf),fill="pink",alpha=0.5) +
  geom_rect(aes(xmin=207.5,xmax=399.5,ymin=-Inf,ymax=Inf),fill="gray",alpha=0.5) +
  geom_line(aes(x=pos,y=cov,group=samp,color=lib)) + theme_classic() + 
  scale_color_manual(values=linecols,labels=linelabs,name="") + 
  xlab("Position (bp)") + ylab("Relative coverage")  + paper.font + 
  theme(legend.key.height=unit(1,"line"),legend.position=c(0.85,0.8)) + ylim(0,0.0053)
dev.off()

# coverage of 3C and genomic libraries overlaid with highlights of region used for ectopic pairing sequence
pdf(paste(out,"/subsequences_coverage_ectopic_paironly.pdf",sep=""),3.3,2)
ggplot(subset(combined.cov,lib!="offtarget"))+ 
  geom_rect(aes(xmin=509.5,xmax=695.5,ymin=-Inf,ymax=Inf),fill="pink",alpha=0.5) +
  geom_line(aes(x=pos,y=cov,group=samp,color=lib)) + theme_classic() + 
  scale_color_manual(values=linecols,labels=linelabs,name="") + 
  xlab("Position (bp)") + ylab("Relative coverage")  + paper.font + 
  theme(legend.key.height=unit(1,"line"),legend.position=c(0.85,0.8)) + ylim(0,0.0053)
dev.off()

# coverage of genomic libraries and expectation
y <- c(c(1:177),rep(178,684),c(177:1))
x <- c(1:1038)
dat <- data.frame(cbind(x,y))
dat$y <- dat$y/sum(dat$y)

pdf(paste(out,"/subsequences_coverage_genomic.pdf",sep=""),3.3,2)
ggplot(subset(combined.cov,lib=="genomic")) + 
  geom_line(aes(x=pos,y=cov,group=samp),color="grey50") + 
  theme_classic() + scale_color_manual(values=linecols,labels=linelabs,name="") + 
  xlab("Position (bp)") + ylab("Relative coverage")  + paper.font + 
  theme(legend.key.height=unit(1,"line"),legend.position=c(0.85,0.8)) + ylim(0,0.0053) +
  geom_line(data=dat,aes(x=x,y=y),linetype="dashed")
dev.off()

