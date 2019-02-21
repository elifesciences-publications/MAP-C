# ScSuRNAseq.R
# plotting and DESeq2 analysis of S. cerevisiae x S. uvarum RNA-seq data
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
library(grid)
library(gplots)
library(scales)
library(DESeq2)

# get gene lengths
genelens <- read.table("nobackup/ScSu_genes_length.txt")
colnames(genelens) <- c("gene", "length")

# samples
samps <- c("ScHAS1pr-TDA1pr_WT_1","ScHAS1pr-TDA1pr_WT_2","ScHAS1pr-TDA1pr_WT_3",
           "ScHAS1pr-TDA1pr_leu3_1","ScHAS1pr-TDA1pr_leu3_2","ScHAS1pr-TDA1pr_leu3_3",
           "ScHAS1pr-TDA1pr_rgt1x2_1","ScHAS1pr-TDA1pr_rgt1x2_2","ScHAS1pr-TDA1pr_rgt1x2_3")

# load data
counts <- genelens
for (samp in samps) {
  temp <- read.table(paste(dir,"/",samp,".counts.txt",sep=""))
  colnames(temp) <- c("gene", samp)
  counts <- merge(counts,temp,by=c("gene"))
}

# extract only raw counts for normalization
rawcounts <- counts[, -(1:2)]
rownames(rawcounts) <- counts$gene

# normalize counts
normed <- sweep(rawcounts, 2, colSums(rawcounts), "/") # by total read count per sample
normed <- sweep(normed, 1, genelens$length, "/") # by gene length (in bp)
normed <- normed*1000000000 # multiply by 1 billion to make FPKM
normwgenes <- cbind(counts$gene, normed) # add back gene annotations
colnames(normwgenes) <- c("gene", samps)

# colors
brewercols <- brewer.pal(5, "Set1")
cols=c("grey",brewercols[c(2, 5)])

# condition names
conds <- c("WT","leu3","rgt1x2")

# function for calculating number of * to add based on p-value
# takes a vector of increasingly stringent (lower) p-value cutoffs and 
# outputs a vector of strings, each with the appropriate number of asterisks
stars <- function(thresh, pval) {
  n = 0
  for (t in thresh) {
    if (pval < t) {
      n = n + 1
    }
  }
  return(paste(rep("*", n), collapse = ""))
}

thresh = c(0.05, 0.01, 0.001, 0.0001) # p-value thresholds for asterisks

# loop through genes
genesofinterest <- c("YMR290C", "YMR291W", "Sbay_13.480", "Sbay_13.481", "YKL038W", "YLR451W")
genelabs <- c("ScHAS1", "ScTDA1", "SuHAS1", "SuTDA1", "ScRGT1", "ScLEU3")
for (i in 1:length(genesofinterest)) {
  geneofinterest <- genesofinterest[i]
  genelab <- genelabs[i]
  
  # subset relevant data
  sub <- subset(normwgenes, gene == geneofinterest)
  table = data.frame
  table <- sub[, 2:4]
  table[2, ] <- sub[, 5:7]
  table[3, ] <- sub[, 8:10]
  colnames(table) <- c("r1", "r2", "r3")
  table$ave <- (table$r1 + table$r2 + table$r3)/3
  table$sd <- apply(table[, 1:3], 1, sd)
  rownames(table) <- conds
  table$cond <- factor(conds, levels = conds)
  
  # calculate p-values 
  gal.test = t.test(table[1, 1:3], table[2, 1:3])
  gal.stars = stars(thresh, gal.test$p.value)
  sat.test = t.test(table[1, 1:3], table[3, 1:3])
  sat.stars = stars(thresh, sat.test$p.value)
  table$stars <- c("", gal.stars, sat.stars)
  
  # bar plot of FPKM of one gene in each condition, with error bars (SEM) and asterisks
  pdf(paste(out,"/ScSuRNAseq_fpkm_", geneofinterest, ".pdf", sep=""), 1.8, 1.5)
  print(ggplot(table) + 
          geom_text(aes(x = cond, y = (ave + sd) + max(table$ave)*.05, label = stars)) + 
          geom_errorbar(aes(x = cond, ymin = ave - sd/sqrt(3), ymax = ave + sd/sqrt(3)), width = .2) + 
          geom_bar(aes(x = cond, fill = cond, y = ave), stat = "identity", color = "black") + 
          theme_classic() + scale_y_continuous(limits = c(0, 1.2*max(table$ave + table$sd)), expand = c(0, 0)) + 
          scale_fill_manual(values = cols) + 
          theme(plot.margin = unit(c(0.2, 0, 0, 0), "cm"), 
                legend.position = "none", 
                plot.title = element_text(hjust = 0.5, face = "italic",size=8), 
                text = element_text(size = 8), axis.text = element_text(size = 8,color="black"),
                axis.text.x = element_text(face="italic")) + 
          xlab("") + ylab("FPKM") + ggtitle(genelab))
  dev.off()
}

# DESeq2 ----

# sample table
samples <- read.table("ScSuRNAseq_table.txt", header=TRUE, stringsAsFactors = TRUE)

# run DESeq2
dds <- DESeqDataSetFromMatrix(countData = rawcounts, 
                              colData = samples, 
                              design = ~ genotype)
dds <- dds[ rowSums(counts(dds)) > 1, ] # filter out genes with <= 1 fragment
dds <- DESeq(dds)

# calculate p-values and log2 Fold Change for leu3 vs. WT
leu <- data.frame(results(dds, contrast = c("genotype", "leu3", "WT")))
leu$gene <- rownames(leu)

# calculate p-values and log2 Fold Change for rgt1x2 vs. WT
rgt <- data.frame(results(dds, contrast = c("genotype", "rgt1x2", "WT")))
rgt$gene <- rownames(rgt)

# volcano plots
#pdf(paste(out,"/ScSuRNAseq_leu3_foldChange_volcano.pdf"), 2.8, 1.6)
png(paste(out,"/ScSuRNAseq_leu3_foldChange_volcano.png",sep=""), 2.8, 1.6, unit="in",res=1200)
ggplot(leu) + geom_point(aes(x = log2FoldChange, y=-log10(pvalue)),size=0.5) + 
  #geom_hline(aes(yintercept=-log10(.05/nrow(genelens))),linetype="dashed") +
  theme_classic() + theme(text = element_text(size = 8), axis.text = element_text(size = 8,color="black"), legend.position = "none") + 
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) + 
  ylab(expression(paste("-",Log[10], " P-value"))) +
  xlab(expression(paste(Log[2], " Fold Change")))
dev.off()

#pdf(paste(out,"/ScSuRNAseq_rgt1x2_foldChange_volcano.pdf"), 2.8, 1.6)
png(paste(out,"/ScSuRNAseq_rgt1x2_foldChange_volcano.png",sep=""), 2.8, 1.6, unit="in",res=1200)
ggplot(rgt) + geom_point(aes(x = log2FoldChange, y=-log10(pvalue)),size=0.5) + 
  #geom_hline(aes(yintercept=-log10(.05/nrow(genelens))),linetype="dashed") +
  geom_text(data=subset(rgt,gene=="YMR291W"),aes(x=log2FoldChange-0.02,y=-log10(pvalue),label="ScTDA1"),hjust=1,size=8*5/14) +
  theme_classic() + theme(text = element_text(size = 8), axis.text = element_text(size = 8,color="black"), legend.position = "none") + 
  scale_y_continuous() + 
  ylab(expression(paste("-",Log[10], " P-value"))) +
  xlab(expression(paste(Log[2], " Fold Change")))
dev.off()

#subset(rgt,pvalue<.05/12510)
