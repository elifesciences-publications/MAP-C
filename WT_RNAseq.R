# WT_RNAseq.R
# plotting of WT RNA-seq data (GSE88952)
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
genelens <- read.table("nobackup/sacCer3_genes_length.txt")
colnames(genelens) <- c("gene", "length")

# samples
samps <- c("asy_r1", "asy_r2", "asy_r3", "gal_r1", "gal_r2", "gal_r3", "sat_r1", "sat_r2", "sat_r3")

# load data
samp.table <- read.table("WT_RNAseq_table.txt", header=T, stringsAsFactors = F)
counts <- read.table(paste(dir,"/asy_r1.counts.txt",sep=""),stringsAsFactors = F)
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

# colors
brewercols <- brewer.pal(4, "Set1")
cols=c("gray",brewercols[1])

# condition names
conds <- c("Exp", "Sat")

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
genesofinterest <- c("YLR451W","YPR022C","YKL038W")
genelabs <- c("LEU3","SDD4","RGT1")

for (i in 1:length(genesofinterest)) {
 geneofinterest <- genesofinterest[i]
 genelab <- genelabs[i]
 
 # subset relevant data
 sub <- subset(normwgenes, gene == geneofinterest)
 table = data.frame
 table <- sub[, 2:4]
 table[2, ] <- sub[, 8:10]
 colnames(table) <- c("r1", "r2", "r3")
 table$ave <- (table$r1 + table$r2 + table$r3)/3
 table$sd <- apply(table[, 1:3], 1, sd)
 rownames(table) <- conds
 table$cond <- factor(conds, levels = conds)

 # calculate p-values 
 sat.test = t.test(table[1, 1:3], table[2, 1:3])
 sat.stars = stars(thresh, sat.test$p.value)
 table$stars <- c("", sat.stars)
 
 # bar plot of FPKM of one gene in each condition, with error bars (SEM) and asterisks
 pdf(paste(out,"/WT_RNAseq_fpkm_", geneofinterest, ".pdf", sep=""), 1.2, 1.2)
 print(ggplot(table) + 
     geom_text(aes(x = cond, y = (ave + sd) + max(table$ave)*.05, label = stars)) + 
     geom_errorbar(aes(x = cond, ymin = ave - sd/sqrt(3), ymax = ave + sd/sqrt(3)), width = .2) + 
     geom_bar(aes(x = cond, fill = cond, y = ave), stat = "identity", color = "black") + 
     theme_classic() + scale_y_continuous(limits = c(0, 1.2*max(table$ave + table$sd/sqrt(3))), expand = c(0, 0)) + 
     scale_fill_manual(values = cols) + 
     theme(plot.margin = unit(c(0.1, 0, 0, 0), "cm"), 
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, face = "italic",size=8), 
        text = element_text(size = 8), axis.text = element_text(size = 8,color="black")) + 
     xlab("") + ylab("FPKM") + ggtitle(genelab))
 dev.off()
}