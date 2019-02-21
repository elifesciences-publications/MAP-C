# ScRNAseq.R
# plotting and DESeq2 analysis of S. cerevisiae RNA-seq data
# Seungsoo Kim

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
library(DESeq2)

# get gene lengths
genelens <- read.table("nobackup/sacCer3_genes_length.txt")
colnames(genelens) <- c("gene", "length")

# load data
samp.table <- read.table("ScRNAseq_table.txt", header=T, stringsAsFactors = F)
counts <- read.table(paste(dir,"/WT_1.counts.txt",sep=""),stringsAsFactors = F)
genes <- counts$V1
counts <- data.frame(counts$V1)
for (samp in samp.table$sample) {
  temp <- read.table(paste(dir,"/",samp,".counts.txt",sep=""))
  colnames(temp) <- c("gene",samp)
  temp <- temp[,2]
  counts <- cbind(counts,temp)
}
colnames(counts) <- c("gene",samp.table$sample)


# exclude non-uniquely assigned read counts
rawcounts <- counts[-((nrow(counts)-4):(nrow(counts))),-1]
rownames(rawcounts) <- counts[1:nrow(rawcounts),]$gene

# normalize counts by total read count per sample
normed <- sweep(rawcounts,2,colSums(rawcounts),"/")
normed <- sweep(normed,1,genelens$length,"/")
normed <- normed*1000000000
normwgenes <- cbind(counts[1:nrow(rawcounts),1],normed)
colnames(normwgenes) <- c("gene",samp.table$sample)

# single gene bar plots ----
# colors
brewercols <- brewer.pal(7,"Set1")
cols=c("grey",brewercols[c(2,3,5)])

# function for calculating number of * to add based on p-value
# takes a vector of increasingly stringent (lower) p-value cutoffs and 
# outputs a vector of strings, each with the appropriate number of asterisks
stars <- function(thresh,pval) {
  n = 0
  for (t in thresh) {
    if (pval < t) {
      n = n + 1
    }
  }
  return(paste(rep("*",n),collapse=""))
}

# loop through genes
genesofinterest <- c("YMR290C","YMR291W","YLR451W","YKL038W","YPR022C")
genelabs <- c("HAS1","TDA1","LEU3","RGT1","SDD4")

for (i in 1:length(genesofinterest)) {
  geneofinterest <- genesofinterest[i]
  genelab <- genelabs[i]
  
  # subset relevant data
  sub <- subset(normwgenes,gene==geneofinterest)
  table = data.frame
  table <- sub[,2:4]
  table[2,] <- sub[,5:7]
  table[3,] <- sub[,8:10]
  table[4,] <- sub[,11:13]
  colnames(table) <- c("r1","r2","r3")
  table$ave <- (table$r1 + table$r2 + table$r3)/3
  table$sd <- apply(table[,1:3],1,sd)
  rownames(table) <- c("WT","leu3D","sdd4D","rgt1D")
  genos <- c("WT","leu3D","sdd4D","rgt1D")
  table$geno <- factor(c("WT","leu3D","sdd4D","rgt1D"),levels=genos)
  
  # calculate p-values 
  thresh = c(.05,.01,.001,.0001)
  leu.test = t.test(table[1,1:3],table[2,1:3])
  leu.stars = stars(thresh,leu.test$p.value)
  sdd.test = t.test(table[1,1:3],table[3,1:3])
  sdd.stars = stars(thresh,sdd.test$p.value)
  rgt.test = t.test(table[1,1:3],table[4,1:3])
  rgt.stars = stars(thresh,rgt.test$p.value)
  table$stars <- c("",leu.stars,sdd.stars,rgt.stars)
  
  pdf(paste(out,"/ScRNAseq_fpkm_",geneofinterest,".pdf",sep=""),2.2,1.8)
  print(ggplot(table) + 
          geom_text(aes(x=geno,y=(ave+sd)+max(table$ave)*.05,label=stars)) + 
          geom_errorbar(aes(x = geno, ymin=ave-sd/sqrt(3), ymax = ave+sd/sqrt(3)),width=.2) + 
          geom_bar(aes(x=geno,fill=geno,y=ave),stat="identity",color="black") + 
          theme_classic() + scale_y_continuous(limits = c(0,1.1*max(table$ave+table$sd)), expand=c(0,0)) + scale_fill_manual(values = cols) + theme(text=element_text(size=8),axis.text=element_text(size=8,color="black"),legend.position="none",plot.title = element_text(hjust=0.5,face ="italic",size=8)) + xlab("") + ylab("FPKM") + ggtitle(genelab))
  dev.off()
}

# DESeq2
samples <- read.table("ScRNAseq_table.txt",header=TRUE,stringsAsFactors = TRUE)

dds <- DESeqDataSetFromMatrix(countData = rawcounts,
                              colData = samples,
                              design = ~ genotype)

dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)

leu <- data.frame(results(dds,contrast=c("genotype","leu3D","WT")))
leu$gene <- rownames(leu)
sdd <- data.frame(results(dds,contrast=c("genotype","sdd4D","WT")))
sdd$gene <- rownames(sdd)
rgt <- data.frame(results(dds,contrast=c("genotype","rgt1D","WT")))
rgt$gene <- rownames(rgt)


gene.names <- read.table("sacCer3_gene_names.txt",header=F,col.names = c("gene","name"))
leu <- merge(leu,gene.names)
sdd <- merge(sdd,gene.names)
sdd[sdd$padj==0 & !is.na(sdd$padj),]$padj <- 10^-38
sdd[sdd$pvalue==0 & !is.na(sdd$pvalue),]$pvalue <- 10^-38
rgt <- merge(rgt,gene.names)
rgt[rgt$padj==0 & !is.na(rgt$padj),]$padj <- 10^-280
rgt[rgt$pvalue==0 & !is.na(rgt$pvalue),]$pvalue <- 10^-280


# volcano plots
pdf(paste(out,"/ScRNAseq_leu3D_foldChange_volcano_padj.pdf",sep=""), 2, 2.5)
ggplot(leu) + geom_point(aes(x=log2FoldChange,y=-log10(padj)),size=0.5) + theme_classic() + 
  geom_text(data=subset(leu,-log10(padj)>20),aes(x=log2FoldChange+0.2,y=-log10(padj),label=name),hjust=0,size=8*5/14) + 
  xlim(min(leu$log2FoldChange),max(leu$log2FoldChange)+1.5) +
  ylab(expression(paste("-", Log[10], " adj. P-value"))) +
  xlab(expression(paste(Log[2], " Fold Change"))) + theme(text=element_text(size=8,color="black"), axis.text=element_text(size=8,color="black"))
dev.off()

pdf(paste(out,"/ScRNAseq_sdd4D_foldChange_volcano_padj.pdf",sep=""), 2, 2.5)
ggplot(sdd) + geom_point(aes(x=log2FoldChange,y=-log10(padj)),size=0.5) + theme_classic() + 
  geom_text(data=subset(sdd,name=="YPR022C"),aes(x=log2FoldChange+0.2,y=-log10(padj),label="SDD4"),hjust=0,size=8*5/14) + 
  geom_text(data=subset(sdd,name=="PCL1"),aes(x=log2FoldChange+0.2,y=-log10(padj)+1,label=name),hjust=0,size=8*5/14) + 
  geom_text(data=subset(sdd,name=="COS8"),aes(x=log2FoldChange+0.05,y=-log10(padj)-1.5,label=name),hjust=0,size=8*5/14) + 
  geom_text(data=subset(sdd,name=="YHL050C"),aes(x=log2FoldChange-0.2,y=-log10(padj),label=name),hjust=1,size=8*5/14) + 
  geom_text(data=subset(sdd,name=="HO"),aes(x=log2FoldChange+0.2,y=-log10(padj),label=name),hjust=0,size=8*5/14) + 
  geom_text(data=subset(sdd,name=="HIS1"),aes(x=log2FoldChange-0.2,y=-log10(padj),label=name),hjust=1,size=8*5/14) + 
  geom_text(data=subset(sdd,name=="RPL33B"),aes(x=log2FoldChange+0.2,y=-log10(padj)+1,label=name),hjust=0,size=8*5/14) + 
  geom_text(data=subset(sdd,name=="ZEO1"),aes(x=log2FoldChange+0.2,y=-log10(padj)-1,label=name),hjust=0,size=8*5/14) + 
  xlim(min(sdd$log2FoldChange),max(sdd$log2FoldChange)+1) +
  ylab(expression(paste("-", Log[10], " adj. P-value"))) +
  xlab(expression(paste(Log[2], " Fold Change"))) + theme(text=element_text(size=8,color="black"), axis.text=element_text(size=8,color="black"))
dev.off()

pdf(paste(out,"/ScRNAseq_rgt1D_foldChange_volcano_padj.pdf",sep=""), 3, 2.5)
ggplot(rgt) + geom_point(aes(x=log2FoldChange,y=-log10(padj)),size=0.5) + theme_classic() + 
  geom_text(data=subset(rgt,-log10(padj)>50 & !(name %in% c("YKR075C","HXT7","HXT4","HXT6","BOP3","IML2","MIG1","YFL054C"))),aes(x=log2FoldChange+0.25,y=-log10(padj),label=name),hjust=0,size=8*5/14) +
  geom_text(data=subset(rgt,(name %in% c("YKR075C","YFL054C"))),aes(x=log2FoldChange-0.25,y=-log10(padj),label=name),hjust=1,size=8*5/14) + 
  geom_text(data=subset(rgt,(name %in% c("BOP3","IML2"))),aes(x=log2FoldChange-0.25,y=-log10(padj)-3,label=name),hjust=1,size=8*5/14) + 
  geom_text(data=subset(rgt,(name %in% c("HXT6"))),aes(x=log2FoldChange-0.2,y=-log10(padj)+8,label=name),hjust=1,size=8*5/14) + 
  geom_text(data=subset(rgt,(name %in% c("MIG1"))),aes(x=log2FoldChange-0.25,y=-log10(padj)+5,label=name),hjust=1,size=8*5/14) + 
  geom_text(data=subset(rgt,(name %in% c("HXT7"))),aes(x=log2FoldChange,y=-log10(padj)+15,label=name),hjust=0.9,size=8*5/14) + 
  geom_text(data=subset(rgt,(name %in% c("HXT4"))),aes(x=log2FoldChange,y=-log10(padj)+15,label=name),hjust=0.1,size=8*5/14) + 
  xlim(min(rgt$log2FoldChange),max(rgt$log2FoldChange)+2) + 
  ylab(expression(paste("-", Log[10], " adj. P-value"))) +
  xlab(expression(paste(Log[2], " Fold Change"))) + theme(text=element_text(size=8,color="black"), axis.text=element_text(size=8,color="black"))
dev.off()

# non adjusted P-value volcano plots
#pdf(paste(out,"/ScRNAseq_leu3D_foldChange_volcano_pvalue.pdf",sep=""), 2, 2.5)
png(paste(out,"/ScRNAseq_leu3D_foldChange_volcano_pvalue.png",sep=""), 2, 2.5, unit="in",res=1200)
ggplot(leu) + geom_point(aes(x=log2FoldChange,y=-log10(pvalue)),size=0.5) + theme_classic() + 
  geom_text(data=subset(leu,-log10(pvalue)>10 & !(name %in% c("ALD5","PCL1"))),aes(x=log2FoldChange+0.2,y=-log10(pvalue),label=name),hjust=0,size=8*5/14) + 
  geom_text(data=subset(leu,-log10(pvalue)>10 & name=="PCL1"),aes(x=log2FoldChange+0.2,y=-log10(pvalue)+5,label=name),hjust=0,size=8*5/14) + 
  geom_text(data=subset(leu,-log10(pvalue)>10 & name=="ALD5"),aes(x=log2FoldChange-0.2,y=-log10(pvalue)-5,label=name),hjust=1,size=8*5/14) + 
  xlim(min(leu$log2FoldChange),max(leu$log2FoldChange)+1.5) +
  ylab(expression(paste("-", Log[10], " P-value"))) +
  xlab(expression(paste(Log[2], " Fold Change"))) + theme(text=element_text(size=8,color="black"), axis.text=element_text(size=8,color="black"))
dev.off()

#pdf(paste(out,"/ScRNAseq_sdd4D_foldChange_volcano_pvalue.pdf",sep=""), 2, 2.5)
png(paste(out,"/ScRNAseq_sdd4D_foldChange_volcano_pvalue.png",sep=""), 2, 2.5, unit="in",res=1200)
ggplot(sdd) + geom_point(aes(x=log2FoldChange,y=-log10(pvalue)),size=0.5) + theme_classic() + 
  geom_text(data=subset(sdd,name=="YPR022C"),aes(x=log2FoldChange+0.2,y=-log10(pvalue),label="SDD4"),hjust=0,size=8*5/14) + 
  geom_text(data=subset(sdd,name=="PCL1"),aes(x=log2FoldChange+0.2,y=-log10(pvalue)+1,label=name),hjust=0,size=8*5/14) + 
  geom_text(data=subset(sdd,name=="COS8"),aes(x=log2FoldChange+0.05,y=-log10(pvalue)-1.5,label=name),hjust=0,size=8*5/14) + 
  geom_text(data=subset(sdd,name=="YHL050C"),aes(x=log2FoldChange-0.2,y=-log10(pvalue),label=name),hjust=1,size=8*5/14) + 
  geom_text(data=subset(sdd,name=="HO"),aes(x=log2FoldChange+0.2,y=-log10(pvalue),label=name),hjust=0,size=8*5/14) + 
  geom_text(data=subset(sdd,name=="HIS1"),aes(x=log2FoldChange-0.2,y=-log10(pvalue),label=name),hjust=1,size=8*5/14) + 
  geom_text(data=subset(sdd,name=="RPL33B"),aes(x=log2FoldChange+0.2,y=-log10(pvalue)+1,label=name),hjust=0,size=8*5/14) + 
  geom_text(data=subset(sdd,name=="ZEO1"),aes(x=log2FoldChange+0.2,y=-log10(pvalue)-1,label=name),hjust=0,size=8*5/14) + 
  xlim(min(sdd$log2FoldChange),max(sdd$log2FoldChange)+1) +
  ylab(expression(paste("-", Log[10], " P-value"))) +
  xlab(expression(paste(Log[2], " Fold Change"))) + theme(text=element_text(size=8,color="black"), axis.text=element_text(size=8,color="black"))
dev.off()

#pdf(paste(out,"/ScRNAseq_rgt1D_foldChange_volcano_pvalue.pdf",sep=""), 3, 2.5)
png(paste(out,"/ScRNAseq_rgt1D_foldChange_volcano_pvalue.png",sep=""), 3, 2.5, unit="in",res=1200)
ggplot(rgt) + geom_point(aes(x=log2FoldChange,y=-log10(pvalue)),size=0.5) + theme_classic() + 
  geom_text(data=subset(rgt,-log10(pvalue)>50 & !(name %in% c("YKR075C","HXT7","HXT4","HXT6","BOP3","IML2","MIG1","YFL054C"))),aes(x=log2FoldChange+0.25,y=-log10(pvalue),label=name),hjust=0,size=8*5/14) +
  geom_text(data=subset(rgt,(name %in% c("YKR075C","YFL054C"))),aes(x=log2FoldChange-0.25,y=-log10(pvalue),label=name),hjust=1,size=8*5/14) + 
  geom_text(data=subset(rgt,(name %in% c("BOP3","IML2"))),aes(x=log2FoldChange-0.25,y=-log10(pvalue)-3,label=name),hjust=1,size=8*5/14) + 
  geom_text(data=subset(rgt,(name %in% c("HXT6"))),aes(x=log2FoldChange-0.2,y=-log10(pvalue)+8,label=name),hjust=1,size=8*5/14) + 
  geom_text(data=subset(rgt,(name %in% c("MIG1"))),aes(x=log2FoldChange-0.25,y=-log10(pvalue)+5,label=name),hjust=1,size=8*5/14) + 
  geom_text(data=subset(rgt,(name %in% c("HXT7"))),aes(x=log2FoldChange,y=-log10(pvalue)+15,label=name),hjust=0.9,size=8*5/14) + 
  geom_text(data=subset(rgt,(name %in% c("HXT4"))),aes(x=log2FoldChange,y=-log10(pvalue)+15,label=name),hjust=0.1,size=8*5/14) + 
  xlim(min(rgt$log2FoldChange),max(rgt$log2FoldChange)+2) + 
  ylab(expression(paste("-", Log[10], " P-value"))) +
  xlab(expression(paste(Log[2], " Fold Change"))) + theme(text=element_text(size=8,color="black"), axis.text=element_text(size=8,color="black"))
dev.off()

# microarray comparisons ----

# load microarray data
leu3ma <- read.table("microarray/leu3D.txt",header=T,stringsAsFactors = F,comment.char = "#",na.strings = "null")
leu3comp <- merge(leu,leu3ma,by.x="gene",by.y="ID_REF")
png(paste(out,"/ScRNAseq_leu3D_sat_v_exp_pvalue.png",sep=""), 2.2, 2.2, units="in",res=1200)
ggplot(leu3comp) + geom_point(aes(x=-log10(P),y=-log10(pvalue),color=name=="LEU3"),size=0.5) + 
  scale_color_manual(values=c("black","red")) + theme_classic() + 
  theme(legend.position="none",text=element_text(size=8,color="black"), axis.text=element_text(size=8,color="black")) +
  ylab(expression(paste("Saturated -", Log[10], " P-value"))) +
  xlab(expression(paste("Exponential -", Log[10], " P-value")))
dev.off()
png(paste(out,"/ScRNAseq_leu3D_sat_v_exp_fc.png",sep=""), 2.2, 2.2, units="in",res=1200)
ggplot(leu3comp) + geom_point(aes(x=VALUE,y=log2FoldChange,color=name=="LEU3"),size=0.5) + 
  scale_color_manual(values=c("black","red")) + theme_classic() + 
  theme(legend.position="none",text=element_text(size=8,color="black"), axis.text=element_text(size=8,color="black")) +
  ylab(expression(paste("Saturated ", Log[2], " Fold Change"))) +
  xlab(expression(paste("Exponential ", Log[2], " Fold Change")))
dev.off()
write.table(x=leu3comp,file = "nobackup/Leu3D_RNAfc.txt")

sdd4ma <- read.table("microarray/sdd4D.txt",header=T,stringsAsFactors = F,comment.char = "#",na.strings = "null")
sdd4comp <- merge(sdd,sdd4ma,by.x="gene",by.y="ID_REF")
png(paste(out,"/ScRNAseq_sdd4D_sat_v_exp_pvalue.png",sep=""), 2.2, 2.2, units="in",res=1200)
ggplot(sdd4comp) + geom_point(aes(x=-log10(P),y=-log10(pvalue),color=name=="YPR022C"),size=0.5) + 
  scale_color_manual(values=c("black","red")) + theme_classic() + 
  theme(legend.position="none",text=element_text(size=8,color="black"), axis.text=element_text(size=8,color="black")) +
  ylab(expression(paste("Saturated -", Log[10], " P-value"))) +
  xlab(expression(paste("Exponential -", Log[10], " P-value")))
dev.off()
png(paste(out,"/ScRNAseq_sdd4D_sat_v_exp_fc.png",sep=""), 2.2, 2.2, units="in",res=1200)
ggplot(sdd4comp) + geom_point(aes(x=VALUE,y=log2FoldChange,color=name=="YPR022C"),size=0.5) + 
  scale_color_manual(values=c("black","red")) + theme_classic() + 
  theme(legend.position="none",text=element_text(size=8,color="black"), axis.text=element_text(size=8,color="black")) +
  ylab(expression(paste("Saturated ", Log[2], " Fold Change"))) +
  xlab(expression(paste("Exponential ", Log[2], " Fold Change")))
dev.off()
write.table(x=sdd4comp,file = "nobackup/Sdd4D_RNAfc.txt")


rgt1ma <- read.table("microarray/rgt1D.txt",header=T,stringsAsFactors = F,comment.char = "#",na.strings = "null")
rgt1comp <- merge(rgt,rgt1ma,by.x="gene",by.y="ID_REF")
png(paste(out,"/ScRNAseq_rgt1D_sat_v_exp_pvalue.png",sep=""), 2.2, 2.2, units="in",res=1200)
ggplot(rgt1comp) + geom_point(aes(x=-log10(P),y=-log10(pvalue),color=name=="RGT1"),size=0.5) + 
  scale_color_manual(values=c("black","red")) + theme_classic() + 
  theme(legend.position="none",text=element_text(size=8,color="black"), axis.text=element_text(size=8,color="black")) +
  ylab(expression(paste("Saturated -", Log[10], " P-value"))) +
  xlab(expression(paste("Exponential -", Log[10], " P-value")))
dev.off()
png(paste(out,"/ScRNAseq_rgt1D_sat_v_exp_fc.png",sep=""), 2.2, 2.2, units="in",res=1200)
ggplot(rgt1comp) + geom_point(aes(x=VALUE,y=log2FoldChange,color=name=="RGT1"),size=0.5) + 
  scale_color_manual(values=c("black","red")) + theme_classic() + 
  theme(legend.position="none",text=element_text(size=8,color="black"), axis.text=element_text(size=8,color="black")) +
  ylab(expression(paste("Saturated ", Log[2], " Fold Change"))) +
  xlab(expression(paste("Exponential ", Log[2], " Fold Change")))
dev.off()
write.table(x=rgt1comp,file = "nobackup/Rgt1D_RNAfc.txt")





rgt1reorder <- rgt[order(rgt$gene=="YMR291W"),]
ggplot(rgt1reorder) + geom_point(aes(x=log2FoldChange,y=-log10(padj),color=gene=="YMR291W"))
rgt1reorder <- rgt[order(rgt$gene=="YMR290C"),]
ggplot(rgt1reorder) + geom_point(aes(x=log2FoldChange,y=-log10(padj),color=gene=="YMR290C"))

sdd4reorder <- sdd[order(sdd$gene=="YMR291W"),]
ggplot(sdd4reorder) + geom_point(aes(x=log2FoldChange,y=-log10(padj),color=gene=="YMR291W"))


subset(normwgenes,gene=="YLR451W")
subset(normwgenes,gene=="YPR022C")
subset(normwgenes,gene=="YKL038W")

subset(leu,padj<.0000001)
subset(sdd,padj<.00000001 & log2FoldChange < -0.9)

rgt1_hits <- subset(rgt,padj<.00000001 & log2FoldChange < -2)







GOI <- "YMR291W"
GOI <- "YBR020W"
counts <- counts[order(counts$gene == GOI),]

ggplot(counts) + geom_point(aes(x=asy_r1,y=asy_r2)) + scale_x_log10() + scale_y_log10()
ggplot(counts) + geom_point(aes(x=asy_r1,y=asy_r3)) + scale_x_log10() + scale_y_log10()
ggplot(counts) + geom_point(aes(x=asy_r1,y=gal_r1,color=gene==GOI)) + scale_x_log10() + scale_y_log10()


genes <- read.table("genes2.info",header=TRUE)
data <- data.frame(genes$gene_id)
samples <- read.table("samples.txt",header=TRUE,stringsAsFactors = FALSE)
for (i in samples$samp) {
  counts <- read.table(paste(i,".counts.txt",sep=""))
  data <- merge(data,counts,by.x="genes.gene_id",by.y="V1")
}
colnames(data) <- c("gene",samples$samp)
countswgenes <- merge(genes,data,by.x="gene_id",by.y="gene")
rownames(data) <- data$gene
data[,1] <- NULL

# joint ----
samples <- read.table("samples2.txt",header=TRUE,stringsAsFactors = FALSE)

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = samples,
                              design = ~ pas + sgn + dg)

dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)

res <- data.frame(results(dds))
res <- data.frame(results(dds, contrast = c("pas","n","y")))
summary(res)

# raw ----

normed <- sweep(countswgenes[,11:34],2,colSums(countswgenes[,11:34]),"/")
normed <- sweep(normed,1,countswgenes$end-countswgenes$start,"/")
normed <- normed*1000000000
normwgenes <- cbind(genes,normed)
