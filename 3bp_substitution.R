# 
# 3bp_substitution.R
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

samps <- c("genomic_1","genomic_2","offtarget_1","offtarget_2","3C_1","3C_2")
brewercols <- brewer.pal(n=6,name="Set1")
paper.font <- theme(text=element_text(size=8), axis.text=element_text(size=8,color="black"), legend.text=element_text(size=8))

# substitution frequency ----
combined <- data.frame(genotype=character())
for (samp in samps) {
  temp <- read.table(paste(dir,"/3bp_substitution_",samp,".subcounts.txt",sep=""),stringsAsFactors=FALSE, col.names=c("genotype","count"))
  temp$samp <- samp

  # normalize by all
  temp$frac <- temp$count/subset(temp,genotype=="All")$count
  temp$count <- NULL

  combined <- rbind(combined,temp)
}
combined <- subset(combined, genotype != "WT" & genotype != "All")
combined$sub <- str_sub(combined$genotype,-1,-1)
combined$pos <- as.numeric(str_sub(combined$genotype,1,str_length(combined$genotype)-1))
casted.total <- cast(combined, pos ~ samp, value="frac", sum)
casted.total$off1 <- casted.total$`offtarget_1`/casted.total$`genomic_1`
casted.total$off2 <- casted.total$`offtarget_2`/casted.total$`genomic_2`
casted.total$off <- (casted.total$off1 + casted.total$off2)/2
casted.total$ratio1 <- casted.total$`3C_1`/casted.total$`genomic_1`
casted.total$ratio2 <- casted.total$`3C_2`/casted.total$`genomic_2`
casted.total$ratio <- (casted.total$ratio1 + casted.total$ratio2)/2

exptstart <- 21
exptend <- 164
exptlen <- 184
offset <- 511

# for paper, Figure 1--figure supplement 2B
pdf(paste(out,"/3bp_substitution_totalperbase.pdf",sep=""),5,2.2)
ggplot(subset(casted.total, pos >= exptstart & pos <= exptend)) + 
  geom_hline(aes(yintercept=0)) +
  geom_segment(aes(x=24.5,xend=34.5,y=1.3,yend=1.3)) +
  geom_segment(aes(x=57.5,xend=64.5,y=1.3,yend=1.3)) +
  geom_segment(aes(x=71.5,xend=81.5,y=1.3,yend=1.3)) +
  geom_segment(aes(x=118.5,xend=128.5,y=1.3,yend=1.3)) +
  geom_segment(aes(x=139.5,xend=149.5,y=1.3,yend=1.3)) +
  geom_segment(aes(x=161.5,xend=168.5,y=1.3,yend=1.3)) +
  geom_text(aes(x=29.5,y=2,label='Leu3'),hjust=0.5,size=2.85) +
  geom_text(aes(x=61,y=2,label='Sdd4'),hjust=0.5,size=2.85) +
  geom_text(aes(x=76.5,y=2,label='Rgt1'),hjust=0.5,size=2.85) +
  geom_text(aes(x=123.5,y=2,label='Rgt1'),hjust=0.5,size=2.85) +
  geom_text(aes(x=144.5,y=2,label='Rgt1'),hjust=0.5,size=2.85) +
  geom_text(aes(x=165,y=2,label='Sdd4'),hjust=0.5,size=2.85) +
  geom_bar(aes(x=pos, y=log2(ratio)), width=1, fill=brewercols[1], stat="identity") + 
  geom_segment(aes(x=pos, xend=pos, y=log2(ratio1), yend=log2(ratio2)), color="black", size=0.2) + 
  theme_classic() + scale_x_continuous(breaks=seq(0,exptlen,25), limits=c(0,exptlen), labels=seq(0,exptlen,25) + offset) + 
  scale_y_continuous(limits=c(-3.4,2.3), breaks=seq(-3,2)) + xlab("Position (bp)") + ylab("Log2 3C/Genomic") + paper.font
dev.off()

# for paper, Figure 1--figure supplement 2C
pdf(paste(out,"/3bp_substitution_totalperbase_offtarget.pdf",sep=""),5,2.2)
ggplot(subset(casted.total, pos >= exptstart & pos <= exptend)) + 
  geom_hline(aes(yintercept=0)) +
  geom_bar(aes(x=pos, y=log2(off)), width=1, fill=brewercols[2], stat="identity") + 
  geom_segment(aes(x=pos, xend=pos, y=log2(off1), yend=log2(off2)), color="black", size=0.2) + 
  theme_classic() + scale_x_continuous(breaks=seq(0,exptlen,25), limits=c(0,exptlen), labels=seq(0,exptlen,25) + offset) + 
  scale_y_continuous(limits=c(-3.4,2.3), breaks=seq(-3,2)) + xlab("Position (bp)") + ylab("Log2 3C/Genomic)") + paper.font
dev.off()

