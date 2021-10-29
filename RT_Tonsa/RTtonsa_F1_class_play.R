#Andrew McCracken
# Copepod Tonsa F1 practice analysis of reciprocal transplant study

## Import or install the libraries that we're likely to need
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)  ### First: BiocManager::install("vsn") AND BiocManager::install("hexbin")

setwd("C:/Users/andre/OneDrive/Desktop/Ecological Genomics/copepods_rnaSeq")

#import the counts matrix
countsTable <- read.table("DE_tonsa_counts_F1.txt", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)
#[1] 24362    16
countsTableRound <- round(countsTable) #Removes Decimals -> becasue DESeq2 doesnt like decimals (and salmons output has decimals)

#import th sample description table
conds <- read.delim("RT_tonsa_F1_samples.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)  
head(conds)

####continued 10/11/21

#lets see how many reads we have from each sample:
colSums(countsTableRound)
mean(colSums(countsTableRound))
barplot(colSums(countsTableRound),names.arg=colnames(countsTableRound),cex.names=0.5,las=3,ytim=c(0,2000000))
abline(h=mean(colSums(countsTableRound)),col="blue",lwd=2) #adds a blue line at the value of the mean number of reads

#the average number of counts per gene
rowSums(countsTableRound) 
mean(rowSums(countsTableRound)) # 11930.81
median(rowSums(countsTableRound)) # 2226

#average number of reads accross all the genes for each of these samples
apply(countsTableRound,2,mean) # 2 in the apply function does the action across columns
apply(countsTableRound,1,mean) # 1 in the apply functon does the action across rows
hist(apply(countsTableRound,1,mean),xlim=c(0,1000),breaks=10000) 

### DESeq --------------------------------------- 

#create a DESeq object and define the experimental design here with a tilda "~"
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, design = ~ line + environment + line:environment)
dim(dds)

# Filter out genes with too few reads - keep reads with average greater than 10 reads/samples
dds <- dds[rowSums(counts(dds)) > 160]  #160 is 10 times our 16 samples
dim(dds)

# Run DESeq model to test for differential expression
dds <- DESeq(dds)

# list the results
resultsNames(dds) 
# "Intercept"   "line_combined_vs_ambient"   "environment_HH_vs_AA"   "linecombined.environmentHH"

#------------------ data Analysis ----------------

### Principle Coordinate Analysis (PCA) to visualize global gene expression patters
vsd <- vst(object = dds, blind = FALSE)
data <- plotPCA(object = vsd, intgroup=c("line", "environment"), returnData=TRUE)
percentVar <- round(100 * attr(data,"percentVar"))

ggplot(data, aes(PC1,PC2,color=environment, shape=line)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

### order and summarize the results from specific contrants
resInteraction <- results(object = dds, alpha = 0.05)
resInteraction <- resInteraction[order(resInteraction$padj),]
head(resInteraction)
tail(resInteraction,1000)

# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
# TRINITY_DN115950_c0_g1   2245.97        4.05237  0.358490   11.3040 1.25396e-29 2.99446e-25
# TRINITY_DN131561_c0_g1   3375.99        4.64570  0.439847   10.5621 4.46620e-26 5.33264e-22
# TRINITY_DN137662_c0_g1  16743.23        4.90200  0.474583   10.3291 5.20658e-25 4.14444e-21
# TRINITY_DN149842_c8_g4  25971.82        4.27274  0.420809   10.1536 3.19275e-24 1.90607e-20
# TRINITY_DN129565_c0_g3  24258.76        4.30553  0.426037   10.1060 5.19661e-24 2.48190e-20
# TRINITY_DN129401_c0_g5  11712.31        4.46355  0.446094   10.0059 1.43650e-23 5.71728e-20

summary(resInteraction)
# out of 24362 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2839, 12%
# LFC < 0 (down)     : 1053, 4.3%
# outliers [1]       : 9, 0.037%
# low counts [2]     : 473, 1.9%
# (mean count < 18)

# About 16% of genes tested show a significant interaction! not differentially expressed genes. for PCA only

resultsNames(resInteraction)


#######################
############################################## TEST FOR EFFECT OF ENVIRONMENT
#######################

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ line + environment)

dds <- DESeq(dds, test="LRT", reduced=~line)
# List the results you've generated

resultsNames(dds)
# [1] "Intercept"  "line_combined_vs_ambient" "environment_HH_vs_AA"


# Order and list and summarize results from specific contrasts
resEnv <- results(dds, alpha = 0.05)
resEnv <- resEnv[order(resEnv$padj),]

head(resEnv)
# DataFrame with 6 rows and 6 columns
# baseMean                    log2FoldChang lfcSE      stat      pvalue            padj
# <numeric>                    <numeric>    <numeric> <numeric>     <numeric>      <numeric>
#   TRINITY_DN138549_c1_g2    582.410       -1.94341  0.173407  118.0825 1.66325e-27 3.96651e-23
# TRINITY_DN138549_c2_g12   773.349       -2.01757  0.203760   91.4210 1.16138e-21 1.38483e-17
# TRINITY_DN150696_c2_g3    297.068        1.31754  0.163636   63.1253 1.93963e-15 1.54188e-11
# TRINITY_DN123676_c0_g2    179.431       -2.51746  0.309813   59.1190 1.48423e-14 7.07917e-11
# TRINITY_DN131329_c1_g1    213.660       -1.23500  0.158361   59.4117 1.27903e-14 7.07917e-11
# TRINITY_DN105043_c0_g1    101.714       -3.94548  0.471847   57.1227 4.09446e-14 1.62741e-10

summary(resEnv)

# 0ut of 24362 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 213, 0.87%
# LFC < 0 (down)     : 235, 0.96%
# outliers [1]       : 41, 0.17%
# low counts [2]     : 473, 1.9%
# (mean count < 18)


resEnv <- resEnv[!is.na(resEnv$padj),]    #get rid of all 'NA' entries excluding empty data
degsEnv <- row.names(resEnv[resEnv$padj < 0.05,]) #order table in order of adjuected p-value
head(resEnv)

#######################
##############################################  TEST FOR EFFECT OF LINE
#######################

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ environment + line)

dds <- DESeq(dds, test="LRT", reduced=~environment)
resultsNames(dds)
# [1] "Intercept"                "environment_HH_vs_AA"     "line_combined_vs_ambient"

resLine <- results(dds, alpha = 0.05)
resLine <- resLine[order(resLine$padj),]

head(resLine)
# log2 fold change (MLE): line combined vs ambient 
# LRT p-value: '~ environment + line' vs '~ environment' 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   TRINITY_DN144155_c0_g8    55.1773       -1.73551  0.199458   75.7313 3.25019e-18 7.13710e-14
# TRINITY_DN116533_c0_g1   101.9489        3.49591  0.435406   53.9339 2.07352e-13 2.27662e-09
# TRINITY_DN142881_c2_g11 1414.7388       -1.43499  0.199146   49.9265 1.59616e-12 1.16834e-08
# TRINITY_DN140379_c0_g5    49.4278        1.69441  0.272190   37.9057 7.42487e-10 4.07607e-06
# TRINITY_DN140379_c0_g6   220.5736        1.86590  0.297107   37.1901 1.07155e-09 4.70602e-06
# TRINITY_DN138009_c0_g2    88.4193        2.12835  0.354530   33.7159 6.37772e-09 2.00069e-05

summary(resLine)
# out of 24362 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 72, 0.3%
# LFC < 0 (down)     : 154, 0.63%
# outliers [1]       : 41, 0.17%
# low counts [2]     : 2362, 9.7%
# (mean count < 25)

resLine <- resLine[!is.na(resLine$padj),]
degsline <- row.names(resLine[resLine$padj < 0.05,])

#######################
##############################################  TEST FOR INTERACTION
#######################

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ environment + line + environment:line)

dds <- DESeq(dds, test="LRT", reduced=~environment + line)
resultsNames(dds)
# [1] "Intercept"                  "environment_HH_vs_AA"       "line_combined_vs_ambient"   "environmentHH.linecombined"


resInt <- results(dds, alpha = 0.05)
resInt <- resInt[order(resInt$padj),]
head(resInt)
# log2 fold change (MLE): environmentHH.linecombined 
# LRT p-value: '~ environment + line + environment:line' vs '~ environment + line' 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   TRINITY_DN115950_c0_g1   2245.97        4.05237  0.358490  118.7642 1.17951e-27 2.76100e-23
# TRINITY_DN131561_c0_g1   3375.99        4.64570  0.439847  101.4199 7.44127e-24 8.70926e-20
# TRINITY_DN137662_c0_g1  16743.23        4.90200  0.474583   95.9225 1.19473e-22 9.32205e-19
# TRINITY_DN149842_c8_g4  25971.82        4.27274  0.420809   94.8761 2.02680e-22 1.18608e-18
# TRINITY_DN129565_c0_g3  24258.76        4.30553  0.426037   93.8792 3.35378e-22 1.57010e-18
# TRINITY_DN129401_c0_g5  11712.31        4.46355  0.446094   91.5220 1.10356e-21 4.30535e-18

summary(resInt)
# out of 24362 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2802, 12%
# LFC < 0 (down)     : 1052, 4.3%
# outliers [1]       : 9, 0.037%
# low counts [2]     : 945, 3.9%
# (mean count < 20)

resInt <- resInt[!is.na(resInt$padj),]
degsInt <- row.names(resInt[resInt$padj < 0.05,])


#######################
##############################################  DATA VISUALIZATION PLOT GENES
#######################

### Plot Individual genes ### 
# Counts of specific top interaction gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds, gene="TRINITY_DN138549_c1_g2", intgroup = (c("line","environment")), returnData=TRUE)
d

p <-ggplot(d, aes(x=environment, y=count, color=line, shape=line, group=line)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p


#######################
############################################## PLOT OVERLAPPING DEGS IN VENN DIAGRAM
#######################
#We'll use the Eulerr package because we all know how nice it is to have the circle scaled. For more info, read here.

library(eulerr)

# Total
length(degsEnv)  # 448
length(degsline)  # 226
length(degsInt)  # 3854

# Intersections
length(intersect(degsEnv,degsline))  # 37
length(intersect(degsEnv,degsInt))  # 44
length(intersect(degsInt,degsline))  # 34

intEL <- intersect(degsEnv,degsline)
length(intersect(degsInt,intEL)) # 7

# Number unique
#Env 448-44-37-7 # 360
#Line 226-37-34-7 # 148
#Interaction 3854-44-34-7 # 3769


fit1 <- euler(c("Env" = 360, "Line" = 148, "Interaction" = 3769, "Env&Line" = 37, "Env&Interaction" = 44, "Line&Interaction" = 34, "Env&Line&Interaction" = 7))

plot(fit1,  lty = 1:3, quantities = TRUE)

plot(fit1, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))

##########
#################### Heatmap of top 20 genes sorted by pvalue
##########

library(pheatmap)

# By interaction

topgenes <- head(rownames(resInt),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df)

# By line

topgenes <- head(rownames(resLine),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df)

#by environment
topgenes <- head(rownames(resEnv),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df)


#Assignment 2 - reciprical environment for 3 generations F3 Data (not tonsa) is there as response to the selection or transplant environment or an interaction between the two. compair the results form class. what does the interation look like? 
#focus on methodological data in deseq
#write out the model of and include the ~

'''
#######################
############################################## TEST FOR EFFECT OF ENVIRONMENT
#######################

dds2 <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ line + environment)
'''


