# HiC.R
# plotting of HiC data
# Seungsoo Kim

setwd("/Volumes/shendure-vol8/projects/mutagenesis.3C/nobackup/MAP-C")

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
#library(Cairo)

# output directory ----
dir.create("figures")
outdir <- "figures"

# color palette
brewercols <- brewer.pal(4,"Set1")

# plot themes ----
paper.fig.rot <- theme(panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank(),
                       panel.background=element_rect(fill="white",color="white"),
                       plot.background=element_rect(fill="white",color="white"),
                       text=element_text(color="black",size=7),
                       axis.text=element_text(color="black",size=7),
                       axis.text.y=element_text(angle=90,hjust=0.5),
                       axis.line=element_line(color="black"),
                       axis.ticks=element_line(color="black"),
                       legend.background=element_rect(fill="white"))

paperhm <- theme(axis.line = element_blank(), 
                 axis.title = element_blank(), 
                 legend.key.size = unit(0.5,"cm"), 
                 legend.title = element_text(size=7), 
                 legend.text = element_text(size=7))


plot_refs <- read.table("plotting_refs.txt",stringsAsFactors=FALSE,sep = "\t")
colnames(plot_refs) <- c("ref","g1","g2","g1p","g2p","rDNA1","rDNA2","bin1","bin2","HAS1TDA1g1","HAS1TDA1g2")

plotting <- read.table("matrices.txt",stringsAsFactors = FALSE)
colnames(plotting) <- c("expt.sample","ref","bsize","mask","filt")
plotting$matr <- paste("data/",plotting$expt.sample,".",plotting$ref,".",plotting$bsize,".matrix.txt",sep="")
plotting$rsum <- paste("data/",plotting$expt.sample,".",plotting$ref,".",plotting$bsize,".rowsums.txt",sep="")

# make figures ----
js <- 1:nrow(plotting)

for (j in js) {
  expt.sample <- plotting[j,]$expt.sample
  print(expt.sample)
  mask <- plotting[j,]$mask
  ref <- plotting[j,]$ref
  bsize <- plotting[j,]$bsize
  filt <- plotting[j,]$filt
  matr <- plotting[j,]$matr
  rsum <- plotting[j,]$rsum
  
  # get reference information
  g1 <- plot_refs[plot_refs$ref==ref,]$g1
  g2 <- plot_refs[plot_refs$ref==ref,]$g2
  g1p <- plot_refs[plot_refs$ref==ref,]$g1p
  g2p <- plot_refs[plot_refs$ref==ref,]$g2p
  nbin1 <- plot_refs[plot_refs$ref==ref,]$bin1
  nbin2 <- plot_refs[plot_refs$ref==ref,]$bin2
  rDNA1 <- plot_refs[plot_refs$ref==ref,]$rDNA1
  rDNA2 <- plot_refs[plot_refs$ref==ref,]$rDNA2
  rDNA <- c(rDNA1,rDNA2)
  
  # set path to save output
  outdir <- "figures"
  
  # load ref annotations
  chrannot <- read.table(paste(ref,".",bsize,".chr_annotations.txt",sep=""))
  binannot <- read.table(paste(ref,".",bsize,".bin_annotations.txt",sep=""))
  
  # load data
  raw <- read.table(matr)
  
  # normalized matrix
  normed.melt <- raw
  if (any(is.na(normed.melt$V4))) {
    normed.melt[is.na(normed.melt$V4),]$V4 <- NA
  }
  
  #add bin annotations
  annotated <- merge(normed.melt,binannot,by.x="V1",by.y="V1")
  colnames(annotated) <- c("bin1","bin2","raw","norm","chr1","cen1","arm1","chrn1")
  annotated <- merge(annotated,binannot,by.x="bin2",by.y="V1")
  colnames(annotated) <- c("bin2","bin1","raw","norm","chr1","cen1","arm1","chrn1","chr2","cen2","arm2","chrn2")
  annotated <- annotated[order(annotated$bin1,annotated$bin2),]
  annotated$bins <- paste(annotated$bin1,annotated$bin2,sep="-")
  annotated$tel1 <- annotated$arm1-annotated$cen1
  annotated$tel2 <- annotated$arm2-annotated$cen2
  
  # heatmap of interactions among HAS1/TDA1 and HXT3
  # for Figure 3--figure supplement 1
  HAS1TDA1g1 <- 292
  HAS1TDA1g2 <- 672
  HXT3g1 <- 80
  HXT3g2 <- 419
  
  subbins <- c(HAS1TDA1g1,HAS1TDA1g2,HXT3g1,HXT3g2)
  annotated.sub <- subset(annotated,bin1 %in% subbins & bin2 %in% subbins)
  annotated.sub$bin1 <- factor(annotated.sub$bin1,levels=subbins)
  annotated.sub$bin2 <- factor(annotated.sub$bin2,levels=subbins)
  
  pdf(paste(outdir,"/HAS1-HXT3_",expt.sample,"_heatmap.pdf",sep=""),1.8,1.8)
  p <- ggplot(annotated.sub) + geom_tile(aes(x=bin1,y=bin2,fill=norm)) +
    paperhm +
    theme_classic() +
    coord_fixed() +
    scale_fill_gradientn(colors=c("white","orange","red","black"), limits=c(0,5), oob=squish, na.value="grey50", name = "O/E") + 
    theme(plot.margin=unit(c(0.2,0.2,0.8,0.8),"cm"),legend.position="none") +
    scale_x_discrete(expand=c(0,0),labels=c("Sc","Su","Sc","Su")) + scale_y_discrete(expand=c(0,0),labels=c("Sc","Su","Sc","Su")) + xlab("") + ylab("") +
    annotation_custom(grob=textGrob(label="HXT3pr",hjust=1,rot=0,gp=gpar(fontsize=8,fontface="italic")),ymin=3.5,ymax=3.5,xmin=-0.6,xmax=-0.6) +
    annotation_custom(grob=textGrob(label="HAS1pr-\nTDA1pr",hjust=1,rot=0,gp=gpar(fontsize=8,fontface="italic")),ymin=1.5,ymax=1.5,xmin=-0.6,xmax=-0.6) +
    annotation_custom(grob=textGrob(label="HXT3pr",hjust=0.5,vjust=1,rot=0,gp=gpar(fontsize=8,fontface="italic")),ymin=-0.4,ymax=-0.4,xmin=3.5,xmax=3.5) +
    annotation_custom(grob=textGrob(label="HAS1pr-\nTDA1pr",hjust=0.5,vjust=1,rot=0,gp=gpar(fontsize=8,fontface="italic")),ymin=-0.4,ymax=-0.4,xmin=1.5,xmax=1.5)
  grid.newpage()
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  chart <- arrangeGrob(gt)
  grid.draw(chart)
  dev.off()
  
  # HXT3-centered heatmap, for Figure 3B
  window <- 5
  subbins1 <- seq(HXT3g1-window,HXT3g1+window)
  subbins2 <- seq(HXT3g2-window,HXT3g2+window)
  annotated.sub <- subset(annotated,bin1 %in% subbins1 & bin2 %in% subbins2)
  
  png(paste(outdir,"/HXT3_",expt.sample,"_heatmap.png",sep=""),1.,1.,"in",res=1200)
  p <- ggplot(annotated.sub) + geom_tile(aes(x=bin1,y=bin2,fill=norm)) + 
    paperhm + 
    coord_fixed() + 
    scale_fill_gradientn(colors=c("white","orange","red","black"), limits=c(0,3), oob=squish, na.value="grey50", name = "O/E") + 
    theme(legend.position="none",text=element_text(size=6)) +
    scale_x_continuous(breaks=c(subbins1[1]-0.5,subbins1[length(subbins1)]+.5),
                       labels=c(bsize/1000*(subbins1[1]-(chrannot[chrannot$V1==paste(g1,"_4",sep=""),]$V2)),bsize/1000*(subbins1[length(subbins1)]-(chrannot[chrannot$V1==paste(g1,"_4",sep=""),]$V2)+1)),
                       expand=c(0,0)) +
    scale_y_continuous(breaks=c(subbins2[1]-0.5,subbins2[length(subbins2)]+.5),
                       labels=c(bsize/1000*(subbins2[1]-(chrannot[chrannot$V1==paste(g2,"_2",sep=""),]$V2)),bsize/1000*(subbins2[length(subbins2)]-(chrannot[chrannot$V1==paste(g2,"_2",sep=""),]$V2)+1)),
                       expand=c(0,0))
  grid.newpage()
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  chart <- arrangeGrob(gt)
  grid.draw(chart)
  grid.force()
  dev.off()
  
  # save comparisons for violin plot
  HXTcomp <- subset(annotated,
                    cen1 >= 15 & cen2 >= 15 & 
                      tel1 > 1 & tel2 > 1 & 
                      bin1 < chrannot[chrannot$V1==paste(g2,"_1",sep=""),2] & 
                      bin2 >= chrannot[chrannot$V1==paste(g2,"_1",sep=""),2] & 
                      bin1 != HAS1TDA1g1 & bin2 != HAS1TDA1g2 &
                      bin1 != HAS1TDA1g1-1 & bin2 != HAS1TDA1g2-1 &
                      bin1 != HAS1TDA1g1+1 & bin2 != HAS1TDA1g2+1 &
                      (chr1 != rDNA1 | chr2 != rDNA2))
  HXTpair <- data.frame(subset(annotated,bin1==HXT3g1 & bin2==HXT3g2)$norm)
  HXTcomp.df <- as.data.frame(HXTcomp$norm)
  
  write.table(HXTcomp.df,paste(outdir,"/HXT3_",expt.sample,"_comp.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep=" ")
  write.table(HXTpair,paste(outdir,"/HXT3_",expt.sample,"_pair.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep=" ")
}

# HXT3 violin plot
all_HXT3 <- data.frame()
filelist <- c(paste(outdir,"/HXT3_ILY456_exponential_Sau3AI",sep=""),
              paste(outdir,"/HXT3_ILY456_exponential_rep2_Sau3AI",sep=""),
              paste(outdir,"/HXT3_YMD3920_exponential_Sau3AI",sep=""),
              paste(outdir,"/HXT3_ILY456_saturated_Sau3AI",sep=""),
              paste(outdir,"/HXT3_YMD3920_saturated_Sau3AI",sep=""),
              paste(outdir,"/HXT3_YMD3266_saturated_Sau3AI",sep=""),
              paste(outdir,"/HXT3_YMD3267_saturated_Sau3AI",sep=""),
              paste(outdir,"/HXT3_YMD3268_saturated_Sau3AI",sep=""),
              paste(outdir,"/HXT3_YMD3269_saturated_Sau3AI",sep="")
)
samplist <- c("exp WT 1",
              "exp WT 2",
              "exp +GATC",
              "sat WT",
              "sat +GATC",
              "sat ymr285-296",
              "sat ymr290",
              "sat ymr290coding",
              "sat ymr290intergenic"
)

for (i in 1:length(filelist)) {
  toadd <- read.table(paste(filelist[i],"_comp.txt",sep=""))
  toadd$V2 <- samplist[i]
  toadd$V3 <- "comp"
  
  toadd_h <- read.table(paste(filelist[i],"_pair.txt",sep=""))
  toadd_h$V2 <- samplist[i]
  toadd_h$V3 <- "pair"
  all_HXT3 <- rbind(all_HXT3,toadd)
  all_HXT3 <- rbind(all_HXT3,toadd_h)
}

# for paper, Figure 3C
pdf(paste(outdir,"/HXT3_violin.pdf",sep=""),2.4,1.4)
ggplot(subset(all_HXT3,V3=="comp")) + 
  geom_violin(aes(x=factor(V2,levels=samplist),y=V1),fill="darkgrey") + 
  paper.fig.rot + 
  geom_point(data=subset(all_HXT3,V3=="pair"),aes(x=factor(V2,levels=samplist),y=V1),color="red",shape="-",size=7)+
  xlab("") + ylab("") + theme(text=element_text(size=8),axis.text.y=element_text(size=8),axis.text.x=element_blank()) +
  scale_x_discrete(labels=samplist)
dev.off()

# labeled
pdf(paste(outdir,"/HXT3_violin_labeled.pdf",sep=""),3.5,3)
ggplot(subset(all_HXT3,V3=="comp")) + 
  geom_violin(aes(x=factor(V2,levels=samplist),y=V1),fill="darkgrey") + 
  paper.fig.rot + 
  geom_point(data=subset(all_HXT3,V3=="pair"),aes(x=factor(V2,levels=samplist),y=V1),color="red",shape="-",size=7)+
  xlab("") + ylab("Normalized interaction frequency") + theme(text=element_text(size=8),axis.text.y=element_text(size=8),axis.text.x=element_text(size=8,angle=90,hjust=0)) +
  scale_x_discrete(labels=samplist)
dev.off()
