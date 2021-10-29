#Andrew McCracken
#Reciprical Transplant Copepod Tonsa F3 Anaysis

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
countsTable_F3 <- read.table("DE_counts_F3.txt", header=TRUE, row.names=1)
head(countsTable_F3)
dim(countsTable_F3) #25279    16
countsTableRound_F3 <- round(countsTable_F3) #Removes Decimals -> becasue DESeq2 doesnt like decimals

#import th sample description table
conds_F3 <- read.delim("RT_tonsa_F3_samples.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1) 
head(conds_F3)


# Reads in each sample:
colSums(countsTableRound_F3)
mean(colSums(countsTableRound_F3))
barplot(colSums(countsTableRound_F3),names.arg=colnames(countsTableRound_F3),cex.names=0.5,las=3,ytim=c(0,2000000))
abline(h=mean(colSums(countsTableRound_F3)),col="blue",lwd=2) #adds a blue line at the value of the mean number of reads

#the average number of counts per gene
rowSums(countsTableRound_F3) 
mean(rowSums(countsTableRound_F3)) # 11220.54
median(rowSums(countsTableRound_F3)) # 2144

#average number of reads accross all the genes for each of these samples
apply(countsTableRound_F3,2,mean) # 2 in the apply function does the action across columns
apply(countsTableRound_F3,1,mean) # 1 in the apply function does the action across rows
hist(apply(countsTableRound_F3,1,mean),xlim=c(0,1000),breaks=10000) 

###--------------------------------DESeq --------------------------------------- 

#create a DESeq object and define the experimental design here with a tilda "~"
dds_F3 <- DESeqDataSetFromMatrix(countData = countsTableRound_F3, colData = conds_F3, design = ~ line + environment + line:environment)
dim(dds_F3) #25279    16

# Filter out genes with too few reads - keep reads with average greater than 10 reads/samples
dds_F3 <- dds_F3[rowSums(counts(dds_F3)) > 160]  #160 is 10 times our 16 samples
dim(dds_F3)

# Run DESeq model to test for differential expression
dds_F3 <- DESeq(dds_F3)

# list the results
resultsNames(dds_F3) 
# "Intercept"   "line_combined_vs_ambient"   "environment_HH_vs_AA"   "linecombined.environmentHH"

#----------------------------data Analysis------------------------------------
vsd_F3 <- vst(object = dds_F3, blind = FALSE)
data_F3 <- plotPCA(object = vsd_F3, intgroup=c("line", "environment"), returnData=TRUE)
percentVar_F3 <- round(100 * attr(data_F3,"percentVar"))

ggplot(data_F3, aes(PC1,PC2,color=environment, shape=line)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ", percentVar_F3[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_F3[2], "% variance")) +
  theme_minimal()

### order and summarize the results from specific contrasts
resInteraction_F3 <- results(object = dds_F3, alpha = 0.05)
resInteraction_F3 <- resInteraction_F3[order(resInteraction_F3$padj),]
head(resInteraction_F3)

# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
# TRINITY_DN142181_c0_g4    384.448        3.11803  0.398309   7.82817 4.95009e-15 1.12976e-10
# TRINITY_DN143012_c0_g4    467.286       -3.31564  0.475434  -6.97391 3.08244e-12 3.51753e-08
# TRINITY_DN131723_c0_g1   1926.655        2.60636  0.387730   6.72209 1.79140e-11 1.36283e-07
# TRINITY_DN142181_c0_g18   364.882        3.11752  0.468239   6.65797 2.77627e-11 1.58407e-07
# TRINITY_DN145818_c5_g1    297.741        1.89854  0.295160   6.43225 1.25730e-10 5.73907e-07
# TRINITY_DN135177_c0_g1   5854.210        2.38301  0.396707   6.00699 1.89001e-09 7.18927e-06

summary(resInteraction_F3)
# out of 25279 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 271, 1.1%
# LFC < 0 (down)     : 60, 0.24%
# outliers [1]       : 5, 0.02%
# low counts [2]     : 2451, 9.7%
# (mean count < 26)

# About 1.3% of genes show a significant interaction -> not differentially expressed genes -> only pulled out one contrast


#######################
############################################## TEST FOR EFFECT OF ENVIRONMENT
#######################

  dds_F3 <- DESeqDataSetFromMatrix(countData = countsTableRound_F3, colData = conds_F3, 
                                 design = ~ line + environment)

dds_F3 <- DESeq(dds_F3, test="LRT", reduced=~line)
# List the results you've generated

resultsNames(dds_F3)
# [1] "Intercept"  "line_combined_vs_ambient" "environment_HH_vs_AA"


# Order and list and summarize results from specific contrasts
resEnv_F3 <- results(dds_F3, alpha = 0.05)
resEnv_F3 <- resEnv_F3[order(resEnv_F3$padj),]

head(resEnv_F3)
# LRT p-value: '~ line + environment' vs '~ line'
# DataFrame with 6 rows and 6 columns
#                           baseMean  log2FoldChange     lfcSE      stat      pvalue        padj
# TRINITY_DN121599_c1_g1     728.155       -6.02233  0.472156  109.3032 1.39269e-25 3.44996e-21
# TRINITY_DN150588_c1_g2    1337.692        2.09424  0.271553   53.2838 2.88676e-13 3.57554e-09
# TRINITY_DN136932_c13_g13   132.738        2.09597  0.291509   48.0239 4.21047e-12 3.47673e-08
# TRINITY_DN146851_c0_g5     465.046        1.33330  0.193383   46.1114 1.11718e-11 6.91867e-08
# TRINITY_DN145745_c0_g5     130.345        1.45951  0.217831   43.6364 3.95424e-11 1.84038e-07
# TRINITY_DN83766_c0_g1      109.739       -1.57117  0.233503   43.4019 4.45756e-11 1.84038e-07

summary(resEnv_F3)
# out of 25279 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 513, 2%
# LFC < 0 (down)     : 315, 1.2%
# outliers [1]       : 16, 0.063%
# low counts [2]     : 491, 1.9%
# (mean count < 19)

# 1.4% of the variance in gene expression comes form the environment


resEnv_F3 <- resEnv_F3[!is.na(resEnv_F3$padj),]    #get rid of all 'NA' entries excluding empty data
degsEnv_F3 <- row.names(resEnv_F3[resEnv_F3$padj < 0.05,]) #order table in order of adjudicated p-value
head(resEnv_F3)

#######################
##############################################  TEST FOR EFFECT OF LINE
#######################

dds_F3 <- DESeqDataSetFromMatrix(countData = countsTableRound_F3, colData = conds_F3, 
                              design = ~ environment + line)

dds_F3 <- DESeq(dds_F3, test="LRT", reduced=~environment)

resultsNames(dds_F3)
# [1] "Intercept"                "environment_HH_vs_AA"     "line_combined_vs_ambient"

resLine_F3 <- results(dds_F3, alpha = 0.05)
resLine_F3 <- resLine_F3[order(resLine_F3$padj),]

head(resLine_F3)
# LRT p-value: '~ environment + line' vs '~ environment' 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
# TRINITY_DN132194_c0_g1  1229.949        1.53901  0.154825   94.6814 2.23631e-22 5.43021e-18
# TRINITY_DN121089_c0_g2   632.855        1.50181  0.163775   80.8863 2.39089e-19 2.90279e-15
# TRINITY_DN134798_c0_g1   152.989       -2.05171  0.224453   77.6336 1.24040e-18 1.00398e-14
# TRINITY_DN129890_c0_g4   153.602       -3.22080  0.341257   76.5916 2.10231e-18 1.27620e-14
# TRINITY_DN147342_c0_g4   151.496        1.86674  0.237816   58.6985 1.83783e-14 8.92522e-11
# TRINITY_DN134960_c1_g9  2869.532        1.94038  0.245758   57.2414 3.85471e-14 1.56000e-10

summary(resLine_F3)
# adjusted p-value < 0.05
# LFC > 0 (up)       : 808, 3.2%
# LFC < 0 (down)     : 837, 3.3%
# outliers [1]       : 16, 0.063%
# low counts [2]     : 981, 3.9%
# (mean count < 21)

#6% of the variance can be explained by line

resLine_F3 <- resLine_F3[!is.na(resLine_F3$padj),]
degsline_F3 <- row.names(resLine_F3[resLine_F3$padj < 0.05,])
head(resLine_F3)

#######################
##############################################  TEST FOR INTERACTION
#######################

dds_F3 <- DESeqDataSetFromMatrix(countData = countsTableRound_F3, colData = conds_F3, 
                              design = ~ environment + line + environment:line)

dds_F3 <- DESeq(dds_F3, test="LRT", reduced=~environment + line)
resultsNames(dds_F3)
# [1] "Intercept"                  "environment_HH_vs_AA"       "line_combined_vs_ambient"   "environmentHH.linecombined"


resInt_F3 <- results(dds_F3, alpha = 0.05)
resInt_F3 <- resInt_F3[order(resInt_F3$padj),]
head(resInt_F3)
# log2 fold change (MLE): environmentHH.linecombined 
# LRT p-value: '~ environment + line + environment:line' vs '~ environment + line' 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
# TRINITY_DN142181_c0_g4    384.448        3.11803  0.398309   59.3373 1.32838e-14 2.96667e-10
# TRINITY_DN143012_c0_g4    467.286       -3.31564  0.475434   46.3684 9.79853e-12 1.09415e-07
# TRINITY_DN131723_c0_g1   1926.655        2.60636  0.387730   43.7887 3.65817e-11 2.72326e-07
# TRINITY_DN142181_c0_g18   364.882        3.11752  0.468239   42.7430 6.24250e-11 3.48535e-07
# TRINITY_DN145818_c5_g1    297.741        1.89854  0.295160   40.8415 1.65090e-10 7.37393e-07
# TRINITY_DN135177_c0_g1   5854.210        2.38301  0.396707   35.1139 3.10977e-09 1.15751e-05

summary(resInt_F3)
# out of 25279 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 235, 0.93%
# LFC < 0 (down)     : 48, 0.19%
# outliers [1]       : 5, 0.02%
# low counts [2]     : 2941, 12%
# (mean count < 28)

# 1% of the variance is explained by interaction

resInt_F3 <- resInt_F3[!is.na(resInt_F3$padj),]
degsInt_F3 <- row.names(resInt_F3[resInt_F3$padj < 0.05,])
head(resInt_F3)

#######################
##############################################  DATA VISUALIZATION PLOT GENES
#######################

### Plot Individual genes ### 
# Counts of specific top interaction gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds_F3, gene="TRINITY_DN132194_c0_g1", intgroup = (c("line","environment")), returnData=TRUE)
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

library(eulerr)

# Total
length(degsEnv_F3)  # 828
length(degsline_F3)  # 1645
length(degsInt_F3)  # 283

# Intersections
length(intersect(degsEnv_F3,degsline_F3))  # 141
length(intersect(degsEnv_F3,degsInt_F3))  # 14
length(intersect(degsInt_F3,degsline_F3))  # 32

intEL_F3 <- intersect(degsEnv_F3,degsline_F3)
length(intersect(degsInt_F3,intEL_F3)) # 7


# Number unique
828-14-141-7 # 666
1645-141-32-7 # 1465
283-14-32-7 # 230


fit1_F3 <- euler(c("Env" = 666, "Line" = 1465, "Interaction" = 230, "Env&Line" = 141, "Env&Interaction" = 14, "Line&Interaction" = 32, "Env&Line&Interaction" = 7))

plot(fit1_F3,  lty = 1:3, quantities = TRUE)

plot(fit1, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))

##########
#################### Heatmap of top 20 genes sorted by pvalue
##########

library(pheatmap)

# By interaction

topgenes_F3 <- head(rownames(resInt_F3),20)
mat_F3 <- assay(vsd_F3)[topgenes_F3,]
mat_F3 <- mat_F3 - rowMeans(mat_F3)
df_F3 <- as.data.frame(colData(dds_F3)[,c("line","environment")])
pheatmap(mat_F3, annotation_col=df_F3)

# By line

topgenes_F3 <- head(rownames(resLine_F3),20)
mat_F3 <- assay(vsd_F3)[topgenes_F3,]
mat_F3 <- mat_F3 - rowMeans(mat_F3)
df_F3 <- as.data.frame(colData(dds_F3)[,c("line","environment")])
pheatmap(mat_F3, annotation_col=df_F3)

#by environment
topgenes_F3 <- head(rownames(resEnv_F3),20)
mat_F3 <- assay(vsd_F3)[topgenes_F3,]
mat_F3 <- mat_F3 - rowMeans(mat_F3)
df_F3 <- as.data.frame(colData(dds_F3)[,c("line","environment")])
pheatmap(mat_F3, annotation_col=df_F3)



























