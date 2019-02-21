# 
# MCAST_hits.R
# 
# Seungsoo Kim
# February 20, 2019

# directories ----
setwd("/Volumes/shendure-vol8/projects/mutagenesis.3C/nobackup/MAP-C")
dir <- "nobackup/mcast"
out <- "figures"

# load libraries
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

# colors
brewercols <- brewer.pal(n=6,name="Set1")
classcols <- c(brewercols[2],brewercols[3],brewercols[5])

# homologous motif clusters with all three motifs in S. cer ----
# for Figure 3
for (pthresh in c("0.001")) {
  for (mgap in c("200")) {
    x <- read.table(paste(dir,"/Scer_",pthresh,"_",mgap,"_Leu3_Sdd4_Rgt1_homologous_LSR.motifs.txt",sep=""))
    x <- x[order(x$V13),]
    x$V12 <- factor(x$V12,levels=c("Leu3","Sdd4","Rgt1"))
    x$genome <- "Scer"
    
    y <- read.table(paste(dir,"/Sbay_",pthresh,"_",mgap,"_Leu3_Sdd4_Rgt1_homologous_LSR.motifs.txt",sep=""))
    y <- y[order(y$V13),]
    y$V12 <- factor(y$V12,levels=c("Leu3","Sdd4","Rgt1"))
    y$genome <- "Suva"
    
    z <- rbind(x,y)
    
    # Figure 3A
    pdf(paste(out,"/MCAST_hits_LSR_",pthresh,"_",mgap,".pdf",sep=""),5,0.5*max(z$V7)+1)
    print(ggplot(z) + geom_rect(aes(xmin=V8,xmax=V9,ymin=0,ymax=-log10(V13),fill=V12),color="black",size=0.3) + 
            geom_segment(aes(x=0,xend=V3-V2,y=0,yend=0)) +
            facet_grid(V7 ~ genome) + theme_classic() + geom_text(aes(x=max(z$V9)/2,y=8,label=V6),hjust=0.5,size=8*5/14,fontface="italic") +
            theme(text=element_text(size=8), axis.text=element_text(size=8,color="black"), legend.text=element_text(size=8,color="black"),
                  panel.spacing = unit(0.35,"cm"), strip.background = element_blank(), 
                  strip.text = element_blank(), 
                  legend.position = c(0.9,1/(max(z$V7)*2)),legend.background = element_blank(), legend.key.size=unit(0.3,"cm")) + 
            xlab("") + ylab("") + scale_fill_manual(values=classcols,name="") + scale_y_continuous(breaks=c(0,3,6),limits=c(0,9)))
    dev.off()
    
    # Figure 3B
    # load ChIP data
    chip <- read.table(paste(dir,"/Scer_",pthresh,"_",mgap,"_Leu3_Sdd4_Rgt1_homologous_LSR_maxFE.txt",sep=""))
    for (s in c("Leu3","Sdd4","Rgt1")) {
      pdf(paste(out,"/MCAST_hits_LSR_ChIP_",s,"_",pthresh,"_",mgap,".pdf",sep=""),0.9,0.5*max(z$V7)+1)
      print(ggplot(subset(chip,V4==s)) + geom_bar(aes(x=V5,y=log2(V3),fill=V5),stat="identity",color="black") + facet_wrap(~V2,ncol=1) + 
              theme_classic() + theme(text=element_text(size=8,color="black"),axis.text=element_text(size=8,color="black"),strip.background = element_blank(),strip.text=element_blank(),legend.position="none") + scale_fill_manual(values=c("grey","red")) + xlab("") + ylab("") + 
              scale_x_discrete(labels=c("Exp","Sat")) + scale_y_continuous(breaks=scales::pretty_breaks(n=3)))
      dev.off()
    }
    
    # Figure 3C
    # motif clusters ranked by total log2 ChIP enrichment
    chip.cast <- cast(subset(chip,V5=="saturated"),V1+V2~V4,value="V3")
    chip.cast <- chip.cast[order(-(log2(chip.cast$Leu3)+log2(chip.cast$Sdd4)+log2(chip.cast$Rgt1))),]
    chip.cast$rank <- 1:nrow(chip.cast)
    pdf(paste(out,"/MCAST_hits_LSR_ChIP_ranked_",pthresh,"_",mgap,".pdf",sep=""),2,2)
    print(ggplot(chip.cast) + geom_bar(aes(x=rank,y=log2(Leu3)+log2(Sdd4)+log2(Rgt1)),stat="identity",color="black") + 
            theme_classic() + theme(text=element_text(size=8,color="black"),
                                    axis.text=element_text(size=8,color="black"),
                                    axis.text.x=element_text(angle=90,vjust=0.5,face="italic"),strip.background = element_blank(),strip.text=element_blank(),legend.position="none") + 
            xlab("") + ylab("") + 
            scale_x_continuous(breaks=chip.cast$rank,labels=chip.cast$V1)+
            scale_y_continuous(breaks=scales::pretty_breaks(n=3)))
    dev.off()
  }
}

# all homologous motif clusters ----
# for Figure 3--figure supplement 1
for (pthresh in c("0.001","0.0005")) {
  for (mgap in c("50","100","200")) {
    x <- read.table(paste(dir,"/Scer_",pthresh,"_",mgap,"_Leu3_Sdd4_Rgt1_homologous.motifs.txt",sep=""))
    x <- x[order(x$V13),]
    x$V12 <- factor(x$V12,levels=c("Leu3","Sdd4","Rgt1"))
    x$genome <- "Scer"
    
    y <- read.table(paste(dir,"/Sbay_",pthresh,"_",mgap,"_Leu3_Sdd4_Rgt1_homologous.motifs.txt",sep=""))
    y <- y[order(y$V13),]
    y$V12 <- factor(y$V12,levels=c("Leu3","Sdd4","Rgt1"))
    y$genome <- "Suva"
    
    z <- rbind(x,y)
    
    # Figure 3--figure supplement 1A
    pdf(paste(out,"/MCAST_hits_",pthresh,"_",mgap,".pdf",sep=""),5,0.5*max(z$V7)+1)
    print(ggplot(z) + geom_rect(aes(xmin=V8,xmax=V9,ymin=0,ymax=-log10(V13),fill=V12),color="black",size=0.3) +
            geom_segment(aes(x=0,xend=V3-V2,y=0,yend=0)) +
            facet_grid(V7 ~ genome) + theme_classic() + geom_text(aes(x=max(z$V9)/2,y=8,label=V6),hjust=0.5,size=7*5/14,fontface="italic") +
            theme(text=element_text(size=8), axis.text=element_text(size=8,color="black"), legend.text=element_text(size=8,color="black"),
                  panel.spacing = unit(0.35,"cm"), strip.background = element_blank(),
                  strip.text = element_blank(),
                  legend.position = c(0.9,1/(max(z$V7)*2)),legend.background = element_blank(), legend.key.size=unit(0.3,"cm")) +
            xlab("") + ylab("") + scale_fill_manual(values=classcols,name="") + scale_y_continuous(breaks=c(0,3,6),limits=c(0,9)))
    dev.off()
    
    # Figure 3--figure supplement 1B
    # load ChIP data
    chip <- read.table(paste(dir,"/Scer_",pthresh,"_",mgap,"_Leu3_Sdd4_Rgt1_homologous_maxFE.txt",sep=""))
    for (s in c("Leu3","Sdd4","Rgt1")) {
      pdf(paste(out,"/MCAST_hits_ChIP_",s,"_",pthresh,"_",mgap,".pdf",sep=""),0.9,0.5*max(z$V7)+1)
      print(ggplot(subset(chip,V4==s)) + geom_bar(aes(x=V5,y=log2(V3),fill=V5),stat="identity",color="black") + facet_wrap(~V2,ncol=1) +
              theme_classic() + theme(text=element_text(size=8,color="black"),axis.text=element_text(size=8,color="black"),strip.background = element_blank(),strip.text=element_blank(),legend.position="none") + scale_fill_manual(values=c("grey","red")) + xlab("") + ylab("") +
              scale_x_discrete(labels=c("Exp","Sat")) + scale_y_continuous(breaks=scales::pretty_breaks(n=3)))
      dev.off()
    }

    # Figure 3--figure supplement 1C
    # motif clusters ranked by total log2 ChIP enrichment
    chip.cast <- cast(subset(chip,V5=="saturated"),V1+V2~V4,value="V3")
    chip.cast <- chip.cast[order(-(log2(chip.cast$Leu3)+log2(chip.cast$Sdd4)+log2(chip.cast$Rgt1))),]
    chip.cast$rank <- 1:nrow(chip.cast)
    pdf(paste(out,"/MCAST_hits_ChIP_ranked_",pthresh,"_",mgap,".pdf",sep=""),3.5,1.7)
    print(ggplot(chip.cast) + geom_bar(aes(x=rank,y=log2(Leu3)+log2(Sdd4)+log2(Rgt1)),stat="identity",color="black") + 
            theme_classic() + theme(text=element_text(size=8,color="black"),axis.text=element_text(size=8,color="black"),axis.text.x=element_text(angle=90,vjust=0.5,face="italic"),strip.background = element_blank(),strip.text=element_blank(),legend.position="none") + 
            xlab("") + ylab("") + 
            scale_x_continuous(breaks=chip.cast$rank,labels=chip.cast$V1)+
            scale_y_continuous(breaks=scales::pretty_breaks(n=3)))
    dev.off()
  }
}

# all motif clusters in each genome ----
for (genome in c("Scer","Sbay")) {
  for (pthresh in c("0.001","0.0005")) {
    for (mgap in c("50","100","200")) {
      x <- read.table(paste(dir,"/",genome,"_",pthresh,"_",mgap,"_Leu3_Sdd4_Rgt1.motifs.txt",sep=""))
      x <- x[order(x$V13),]
      x$V12 <- factor(x$V12,levels=c("Leu3","Sdd4","Rgt1"))
      
      pdf(paste(out,"/MCAST_hits_",genome,"_",pthresh,"_",mgap,".pdf",sep=""),2.5,0.5*max(x$V7)+1)
      print(ggplot(x) + geom_rect(aes(xmin=V8,xmax=V9,ymin=0,ymax=-log10(V13),fill=V12),color="black",size=0.3) + 
              geom_segment(aes(x=0,xend=V3-V2,y=0,yend=0)) +
              facet_wrap(~V7,ncol=1) + theme_classic() + geom_text(aes(x=max(x$V9)/2,y=8,label=V6),hjust=0.5,size=7*5/14) +
              theme(text=element_text(size=8), axis.text=element_text(size=8,color="black"), legend.text=element_text(size=8,color="black"),
                    panel.spacing = unit(0.35,"cm"), strip.background = element_blank(), 
                    strip.text = element_blank(), 
                    legend.position = c(0.9,0.04),legend.background = element_blank(), legend.key.size=unit(0.3,"cm")) + 
              xlab("") + ylab("") + scale_fill_manual(values=classcols,name="") + scale_y_continuous(breaks=c(0,3,6),limits=c(0,9)))
      dev.off()
    }
  }
}