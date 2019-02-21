# 
# in-gene_ko.R
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

# for paper, Figure 2--figure supplement 1A ----
height = .3

gnames <- read.table("references/sacCer3_gene_names.txt")
cen <- read.table("references/centromere_positions.txt")
tel <- read.table("references/telomere_positions.txt")
genes <- read.table("in-gene_ko_TF_positions.txt")
ctlgenes <- read.table("in-gene_ko_positions.txt")
chr <- read.table("references/chr_positions.txt")
cen$V5 <- "+"
tel$V5 <- "+"
chr$V5 <- "+"

cen$V6 <- "centromere"
tel$V6 <- "telomere"
ctlgenes$V6 <- "TFKO ctl gene"
genes$V6 <- "TFKO gene"
chr$V6 <- "chromosome"

chroms <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI","MT")
chroms <- paste("chr",chroms,sep="")

cen$V2 <- factor(cen$V2,levels=chroms)
tel$V2 <- factor(tel$V2,levels=chroms)
ctlgenes$V2 <- factor(ctlgenes$V2,levels=chroms)
genes$V2 <- factor(genes$V2,levels=chroms)
chr$V2 <- factor(chr$V2,levels=chroms)

combined <- rbind(cen,tel,ctlgenes,genes,chr)
combined$V6 <- factor(combined$V6,levels=c("telomere","TFKO gene","centromere","chromosome","TFKO ctl gene"))


ovalFun <- function(center = c(0,0),height = 1, width = 1, npoints = 100){
  r = height / 2
  r2 = width / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r2 * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

circ1 <- ovalFun(c((cen[1,3]+cen[1,4])/2,16-1+height/2),height,20000,100)
circ2 <- ovalFun(c((cen[2,3]+cen[2,4])/2,16-2+height/2),height,20000,100)
circ3 <- ovalFun(c((cen[3,3]+cen[3,4])/2,16-3+height/2),height,20000,100)
circ4 <- ovalFun(c((cen[4,3]+cen[4,4])/2,16-4+height/2),height,20000,100)
circ5 <- ovalFun(c((cen[5,3]+cen[5,4])/2,16-5+height/2),height,20000,100)
circ6 <- ovalFun(c((cen[6,3]+cen[6,4])/2,16-6+height/2),height,20000,100)
circ7 <- ovalFun(c((cen[7,3]+cen[7,4])/2,16-7+height/2),height,20000,100)
circ8 <- ovalFun(c((cen[8,3]+cen[8,4])/2,16-8+height/2),height,20000,100)
circ9 <- ovalFun(c((cen[9,3]+cen[9,4])/2,16-9+height/2),height,20000,100)
circ10 <- ovalFun(c((cen[10,3]+cen[10,4])/2,16-10+height/2),height,20000,100)
circ11 <- ovalFun(c((cen[11,3]+cen[11,4])/2,16-11+height/2),height,20000,100)
circ12 <- ovalFun(c((cen[12,3]+cen[12,4])/2,16-12+height/2),height,20000,100)
circ13 <- ovalFun(c((cen[13,3]+cen[13,4])/2,16-13+height/2),height,20000,100)
circ14 <- ovalFun(c((cen[14,3]+cen[14,4])/2,16-14+height/2),height,20000,100)
circ15 <- ovalFun(c((cen[15,3]+cen[15,4])/2,16-15+height/2),height,20000,100)
circ16 <- ovalFun(c((cen[16,3]+cen[16,4])/2,16-16+height/2),height,20000,100)

tfs <- merge(gnames,subset(combined,V6=="TFKO gene"),all.y=TRUE,by.x="V1",by.y="V1")
colnames(tfs) <- c("id","gene","chr","start","end","dir","class")

pdf(paste(out,"/in-gene_ko_schematic.pdf",sep=""),4.2,3.8)
ggplot(combined) + 
  geom_text(aes(x=852110,y=3+height/2-.1,label="*"),size=5*8/14) +
  geom_text(aes(x=852110,y=3+height+.3,label="HAS1pr-TDA1pr"),size=5*8/14,fontface="italic") +
  geom_text(data=subset(tfs,!(gene %in% c("MSN4","RGT1"))),aes(x=(start+end)/2,y=16-as.numeric(chr)+height+.3,label=gene),size=5*8/14,fontface="italic") +
  geom_text(data=subset(tfs,(gene %in% c("MSN4"))),aes(x=(start+end)/2-60000,y=16-as.numeric(chr)+height+.3,label=gene),size=5*8/14,fontface="italic") +
  geom_text(data=subset(tfs,(gene %in% c("RGT1"))),aes(x=(start+end)/2+60000,y=16-as.numeric(chr)+height+.3,label=gene),size=5*8/14,fontface="italic") +
  geom_rect(data=subset(combined, V6=="telomere"), aes(xmax=V4, xmin=V3, ymax=16-as.numeric(V2)+height, ymin=16-as.numeric(V2)), fill="gray") +
  geom_polygon(data=circ1,aes(x,y),fill="gray") +
  geom_polygon(data=circ2,aes(x,y),fill="gray") +  
  geom_polygon(data=circ3,aes(x,y),fill="gray") +  
  geom_polygon(data=circ4,aes(x,y),fill="gray") +  
  geom_polygon(data=circ5,aes(x,y),fill="gray") +
  geom_polygon(data=circ6,aes(x,y),fill="gray") +  
  geom_polygon(data=circ7,aes(x,y),fill="gray") +  
  geom_polygon(data=circ8,aes(x,y),fill="gray") +  
  geom_polygon(data=circ9,aes(x,y),fill="gray") +  
  geom_polygon(data=circ10,aes(x,y),fill="gray") +  
  geom_polygon(data=circ11,aes(x,y),fill="gray") +  
  geom_polygon(data=circ12,aes(x,y),fill="gray") +  
  geom_polygon(data=circ13,aes(x,y),fill="gray") +  
  geom_polygon(data=circ14,aes(x,y),fill="gray") +  
  geom_polygon(data=circ15,aes(x,y),fill="gray") +  
  geom_polygon(data=circ16,aes(x,y),fill="gray") +  
  geom_segment(data=subset(combined, V6=="TFKO ctl gene"),aes(x=(V3+V4)/2, xend=(V3+V4)/2, y=16-as.numeric(V2)+height, yend=16-as.numeric(V2)),color="black",size=0.2) + 
  geom_segment(data=subset(combined, V6=="TFKO gene"),aes(x=(V3+V4)/2, xend=(V3+V4)/2, y=16-as.numeric(V2)+height, yend=16-as.numeric(V2)),color="red",size=0.2) + 
  geom_rect(data=subset(combined, V6=="chromosome"), aes(xmax=V4, xmin=V3, ymax=16-as.numeric(V2)+height, ymin=16-as.numeric(V2)),color="black",fill="white",alpha=0,size=0.2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line.x = element_line(color="black"), axis.line.y = element_blank(), 
                     axis.text.x = element_text(size=8,color="black"), axis.text.y = element_text(size=8, vjust=0.5, hjust=0, color="black"), 
                     axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_text(size=8,color="black")) + 
  xlab("Genomic coordinate (kb)") + scale_y_discrete(limits=seq(0,15)+height/2,labels=rev(chroms[-17])) + 
  scale_x_continuous(breaks=c(0,500000,1000000,1500000),labels=c(0,500,1000,1500))
dev.off()

# ----
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

# for paper, Figure 2--figure supplement 1B

annot <- read.table("in-gene_ko_annotations.txt")
m <- merge(groups,annot,by.x="gene",by.y="V5")
m <- merge(m,gnames,by.x="gene",by.y="V1")
colnames(m) <- c("gene","group","chr","start","end","strand","name")
m$namepos <- 1
m[m$name %in% c("EBS1","YDR209C","COQ4","GAL83","CHZ1","YPL247C","YPL245W","YOR111W","YOR114W","SUB1","FPR4","RIF2","YKL063C","AGA2","YGL034C","NUP120","BLI1","NFU1","VAM7","SIP2","YPT32","YKL033W-A"),]$namepos <- 2
m$adjust=0
m[m$name == "VAM7",]$adjust=-1
for (g in unique(m$group)) {
  print(g)
  x <- subset(m,group == g)
  x$dir <- as.numeric(x$strand=="+")
  pdf(paste(out,"/in-gene_ko_group_",g,".pdf",sep=""),3.2,0.8)
  print(ggplot(x) + geom_rect(aes(xmin=start/1000,xmax=end/1000,ymin=0,ymax=0.5,fill=gene==g)) + 
          geom_text(aes(x=(start+end)/2000+adjust,y=namepos,label=name),size=7*5/14, fontface="italic") +
          scale_fill_manual(values=c("black","red")) + theme_classic() + paper.font +
          theme(axis.line.y = element_blank(),axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position="none") +
          xlab("Chromosome coordinate (kb)") + ylab("") + ylim(0,2.5) + scale_x_continuous(expand=c(0.05,0.05)))
  dev.off()
}


# for paper, Figure 2--figure supplement 1C
combined$groupfac <- factor(combined$group)
levels(combined$groupfac) <- c("UME6","MIG3","MIG1","MIG2","RGT1","MSN4","LEU3","MSN2","AZF1","GAL4")

pdf(paste(out,"/in-gene_ko_effect_of_centromere_distance.pdf",sep=""),2.7,4)
ggplot(subset(combined,`genomic_1` > .002 & `genomic_2` > .002)) + 
  geom_point(aes(x=cendist/1000,y=ratio1,color=gene!=group)) + 
  geom_point(aes(x=cendist/1000,y=ratio2,color=gene!=group)) + 
  geom_segment(aes(x=cendist/1000,xend=cendist/1000,y=ratio1,yend=ratio2,color=gene!=group)) + 
  theme_classic() + scale_color_manual(values=c("red","black"),labels=c("TF","Neighboring gene"),name="") + 
  xlab("Distance from centromere (kb)") + ylab("3C/Genomic") + 
  theme(text=element_text(size=8,color="black"),legend.position=c(0.75,0.1),legend.background = element_blank()) + paper.font
dev.off()


# for paper, Figure 2--figure supplement 1D
pdf(paste(out,"/in-gene_ko_pairing_ratio_v_input_bygroup.pdf",sep=""),6.6,3.5)
ggplot(combined[order(combined$gene!=combined$group),]) + 
  geom_point(aes(x=log10(`genomic_1`),y=log2(ratio1),color=gene!=group)) + 
  geom_point(aes(x=log10(`genomic_2`),y=log2(ratio2),color=gene!=group)) + 
  geom_segment(aes(x=log10(`genomic_1`),xend=log10(`genomic_2`),y=log2(ratio1),yend=log2(ratio2),color=gene!=group)) + 
  theme_bw() + scale_x_continuous(limits = c(-5.5,-.5)) + 
  facet_wrap(~groupfac,ncol=5) + 
  scale_color_manual(values=c("red","black"),labels=c("TF","Neighboring gene"),name="") + 
  ylab(expression(paste(Log[2], " 3C/Genomic"))) + xlab(expression(paste(Log[10]," Frequency in genomic"))) + 
  scale_y_continuous() + theme(strip.background=element_blank(),strip.text=element_text(size=8,face="italic")) + paper.font
dev.off()

