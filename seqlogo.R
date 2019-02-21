# 
# seqlogo.R
# 
# Seungsoo Kim
# February 20, 2019

# directories ----
setwd("/Volumes/shendure-vol8/projects/mutagenesis.3C/nobackup/MAP-C")
dir <- "yetfasco"
out <- "figures"

# load libraries ----
#install.packages("ggseqlogo")
library("ggplot2")
library("ggseqlogo")

# list of motif names
motifs <- read.table("yetfasco.txt")

for (i in motifs$V1) {
  # forward motif
  # read PFM
  x <- read.table(paste(dir,"/",i,"_pad15.pfm",sep=""))
  rownames(x) <- x$V1
  x$V1 <- NULL
  x <- as.matrix(x)
  # make sequence logo
  pdf(paste(out,"/seqlogo_",i,".pdf",sep=""),1.8,0.7)
  print(ggseqlogo(x) + theme(axis.text.x=element_blank()) + scale_y_continuous(breaks=c(0,1,2),limits=c(0,2)) + theme(axis.line.y = element_line(), axis.ticks.y = element_line(), text=element_text(size=8)))
  dev.off()
  
  # reverse complement
  # read PFM
  x <- read.table(paste(dir,"/",i,"_rc_pad15.pfm",sep=""))
  rownames(x) <- x$V1
  x$V1 <- NULL
  x <- as.matrix(x)

  # make sequence logo
  pdf(paste(out,"/seqlogo_",i,"_rc.pdf",sep=""),1.8,0.7)
  print(ggseqlogo(x) + theme(axis.text.x=element_blank()) + scale_y_continuous(breaks=c(0,1,2),limits=c(0,2)) + theme(axis.line.y = element_line(), axis.ticks.y = element_line(), text=element_text(size=8)))
  dev.off()
}
