#hwk 3  
library(ggplot2)
library(gridExtra)
library(ggbiplot)


setwd("C:/Users/andre/OneDrive/Desktop/Ecological Genomics/Pop-Land-Genomics")
# Get the list of admixed individuals:
Admixed <- read.table("Admixed.Inds",header=F)

# Get the meta data:
meta <- read.table("Combined_Transect_Sampling_Data_2020.txt", sep="\t",header=T)

# merge them together:
meta_admx <- merge(meta, Admixed, by.x="ID", by.y="V1")
str(meta_admx)  

# Read in the Admixture coefficients for KBals that we made from the K=5 file:
KBals <- read.table("Admixed_KBals", sep="\t", header=F)
names(KBals) = c("ID","KBals")

# Second merge:
meta_admx_KBals <- merge(meta_admx,KBals,by="ID")


# Bring in phenotype data:
pheno <- read.table("VT_Garden_Phenotypes_2021.txt",sep="\t",header=T)
clim <- read.table("climDat.txt",sep="\t",header=T)
# Merge pheno data with meta and KBals:
meta_admx_KBals_clim <- merge(meta_admx_KBals,clim,by="ID")
meta_admx_KBals_cp <- merge(meta_admx_KBals_clim,pheno,by="ID")
# This is the average date of the last freezing event in spring, after which temperatures stay above 0C for the rest of the growing season.
plotmean <- ggplot(meta_admx_KBals_cp,aes(x=KBals,y=mean_finalFreeze, color=Transect.x, ellipse = TRUE)) +
  geom_point(size=2) +
  xlab("Proportion P. balsamifera ancestry") +
  ylab("final freeze") 

plotmean

# This is the average number of growing degree days (a measure of spring warming) that accumulate from Jan01 to the date of last freeze in spring.
plotgdd <- ggplot(meta_admx_KBals_cp,aes(x=KBals,y=mean_cGDDfreeze, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("mean_cGDDfreeze") 

plotgdd

#This is the mean number of chilling degree days across the year.  Chilling degree days are thought to be a cue that plants use to determine how much winter they have passed through, which can be used to decide when to break dormancy (flush).
plotmed_DD0 <- ggplot(meta_admx_KBals_cp,aes(x=KBals,y=med_DD0, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("med_DD0") 

#individuals with trichocarp ancestry take longer in spring to break bud, higher balsamifera== earlier
#populations that have evolved under warmer climates require more heat to bud
#natural pops in warm know that they can get a few warm days even in winter

plotmed_DD0

grid.arrange(plotmed_DD0, plotmean, plotgdd, nrow = 3)


#bud set for trait
budflush <- read.table("plink2.FLUSH.glm.linear",skip=1,sep="\t",header=F)
names(budflush) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
budflush <- budflush[which(budflush$TEST=="ADD"),]



# regression plots
snps <- read.table("Chr03.kept.sites", sep="\t", heade=T)

budflush2 <- cbind(snps, budflush[,-c(1,2)])
budflush2$outlier = ifelse(budflush2$P<quantile(budflush2$P,0.001),2,1)
#regressions for each of these plots

# linear models testing trait ~ genome-wide admixture association
#mean final freeze
summary(lm(mean_finalFreeze~KBals + Transect.x, data=meta_admx_KBals_cp))

summary(lm(mean_cGDDfreeze~KBals + Transect.x, data=meta_admx_KBals_cp))

summary(lm(med_DD0~KBals + Transect.x, data=meta_admx_KBals_cp))

lm(formula = med_DD0 ~ KBals + Transect.x, data = meta_admx_KBals_cp)

# What about the effects of local ancestry within the genome, after controlling for genome-wide ancestry effects as a covariate in the GLM model?



######  Bring in Association results from Plink   ######


# Define association outliers as the upper 0.1% of p-values
#define snps

####################
budflush <- read.table("plink2.FLUSH.glm.linear",skip=1,sep="\t",header=F)
names(budflush) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
budflush <- budflush[which(budflush$TEST=="ADD"),]
budflush2 <- cbind(snps, budflush[,-c(1,2)])
budflush2$outlier = ifelse(budflush2$P<quantile(budflush2$P,0.01),2,1)

pflush <- ggplot(budflush2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=budflush2$outlier, color=budflush2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Bud flush")

#pflush




ffreeze <- read.table("plink2.mean_finalFreeze.glm.linear",skip=1,sep="\t",header=F)
names(ffreeze) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
ffreeze2 <- ffreeze[which(ffreeze$TEST=="ADD"),]
head(ffreeze2)

#########  ffreeze2  #########

ffreeze2 <- cbind(snps, ffreeze2[,-c(1:2)])
ffreeze2$outlier = ifelse(ffreeze2$P<quantile(ffreeze2$P,0.01),2,1) #upper .1% 

p1 <- ggplot(ffreeze2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=ffreeze2$outlier, color=ffreeze2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("final freeze")

#p1

####### GDD #########
gdd <- read.table("plink2.mean_cGDDfreeze.glm.linear",skip=1,sep="\t",header=F)
names(gdd) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
gdd <- gdd[which(gdd$TEST=="ADD"),]
gdd2 <- cbind(snps, gdd[,-c(1,2)])
gdd2$outlier = ifelse(gdd2$P<quantile(gdd2$P,0.01),2,1)

p2 <- ggplot(gdd2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=gdd2$outlier, color=gdd2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("gdd")

#p2

#########  med DD0  #########
dd0 <- read.table("plink2.med_DD0.glm.linear",skip=1,sep="\t",header=F)
names(dd0) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
dd0 <- dd0[which(dd0$TEST=="ADD"),]
dd02 <- cbind(snps, dd0[,-c(1,2)])
dd02$outlier = ifelse(dd02$P<quantile(dd02$P,0.01),2,1)

p3 <- ggplot(dd02,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=dd02$outlier, color=dd02$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("dd0")

#p3

#grid.arrange(p1, p2, p3, pflush, nrow = 4)

# -------------

# Get outliers for a given trait association:
ffreeze_outliers <- ffreeze2[which(ffreeze2$outlier==2),c(2,3,9)]
gdd_outliers <- gdd2[which(gdd2$outlier==2),c(2,3,9)]
dd0_outliers <- dd02[which(dd02$outlier==2),c(2,3,9)]
budflush_outliers <- budflush2[which(budflush2$outlier==2),c(2,3,9)]

#plot ancestry
# Read in list of positions
snps <- read.table("Chr03.kept.sites",sep="\t", header=T)

# Plot freq of LAI along chr
AF <- read.table("Chr03_LAI_freq.afreq", skip=1,sep="\t",header=F)
names(AF) = c("CHROM",  "ID",   "REF",  "ALT",  "ALT_FREQS",    "OBS_CT")
str(AF)

AF2 <- cbind(snps,AF)

windows <- seq(1,max(AF2$POS),5e4)
AF_windows <- numeric()

for(i in 1:length(windows)){
  tmp=AF2[which(AF2$POS>windows[i] & AF2$POS<windows[i+1]),"ALT_FREQS"]
  ancfreq=mean(tmp)
  AF_windows[i] = ancfreq
}

AF3 <- as.data.frame(cbind(windows,AF_windows))
names(AF3) = c("window","AvgAncFreq")
#
upper = mean(AF3$AvgAncFreq,na.rm=T) + 2*sd(AF3$AvgAncFreq,na.rm=T)
lower = mean(AF3$AvgAncFreq,na.rm=T) - 2*sd(AF3$AvgAncFreq,na.rm=T)

outliers_upper = AF3[which(AF3$AvgAncFreq>upper),]
outliers_lower = AF3[which(AF3$AvgAncFreq<lower),]

# Print the outlier regions out
outliers_upper
outliers_lower

# And finally, make the 4-panel plot with the trait associations
p4 <- ggplot(AF3[,-3],aes(x=window,y=AvgAncFreq)) +
  geom_line(size=0.8, color="blue") + 
  xlab("Position (bp) along chromosome") +
  ylab("Frequency P. trichocarpa ancestry") +
  geom_hline(yintercept=mean(AF2$ALT_FREQS), color = "red") + 
  geom_hline(yintercept=upper, linetype="dashed", color = "red") + 
  geom_hline(yintercept=lower, linetype="dashed", color = "red") +
  ggtitle("Chr03: Local ancestry")

#p4

#grid.arrange(p1, p2, p3, pflush, p4, nrow = 5)
#grid.arrange(p1, p2, p3, p4, nrow = 4)

# Get the betas from each trait and look at pleiotropy between traits
betas <- cbind(ffreeze2[,c(1:3,9)],gdd2[,9],dd02[,9],budflush2[,9])
names(betas) = c("CHROM","POS","ID","beta_ffreeze","beta_gdd","beta_dd0", "beta_flush")
str(betas)

cor(betas[,4:6],betas[4:6])

plot(beta$beta_ffreeze,betas$beta_flush)

p5 <- ggplot(betas,aes(x=beta_ffreeze,y=beta_flush)) +
  geom_point(color="darkgray") + 
  xlab("first freeze") +
  ylab("bud flush") +
  ggtitle("Correlation of first freeze and bud flush")

#p5

p6 <- ggplot(betas,aes(x=beta_gdd,y=beta_flush)) +
  geom_point(color="darkgray") + 
  xlab("gdd") +
  ylab("bud flush") +
  ggtitle("Correlation of gdd and bud flush")

#p6

p7 <- ggplot(betas,aes(x=beta_dd0,y=beta_flush)) +
  geom_point(color="darkgray") + 
  xlab("dd0") +
  ylab("bud flush") +
  ggtitle("Correlation of dd0 and bud flush")

#p7

grid.arrange(p5, p6, p7, nrow = 3)

#__________________outliers overlaps______________________

library(GenomicRanges)

# Calculate Genomic Ranges for outlier bins
CHR="Chr03"
# Define the genomic ranges of the budflush bins
budflushGR <- GRanges(CHR,IRanges(budflush2$POS-2.5e4,budflush2$POS+2.5e4),POS=budflush2$POS, P=budflush2$P, outlier=budflush2$outlier)

budflushGRout <- unlist(reduce(split(budflushGR, ~outlier)))
budflushGRout$outlier <- names(budflushGRout)
budflushGRCand <- subset(budflushGRout, outlier==2)

budflushGRCand # Print the candidate regions

# GRanges object with 6 ranges and 1 metadata column:
#   seqnames            ranges strand |     outlier
# <Rle>         <IRanges>  <Rle> | <character>
#   2    Chr03   9031506-9100590      * |           2
# 2    Chr03 11993834-12154546      * |           2
# 2    Chr03 12447851-12523308      * |           2
# 2    Chr03 18515267-18776502      * |           2
# 2    Chr03 18778201-18861019      * |           2
# 2    Chr03 19698322-19768188      * |           2
# -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths


####### Define the genomic ranges of the ddo bins
dd0GR <- GRanges(CHR,IRanges(dd02$POS-2.5e4,dd02$POS+2.5e4),POS=dd02$POS, P=dd02$P, outlier=dd02$outlier)

dd0GRout <- unlist(reduce(split(dd0GR, ~outlier)))
dd0GRout$outlier <- names(dd0GRout)
dd0GRCand <- subset(dd0GRout, outlier==2)

dd0GRCand # Print the candidate regions

# GRanges object with 8 ranges and 1 metadata column:
#   seqnames            ranges strand |     outlier
# <Rle>         <IRanges>  <Rle> | <character>
#   2    Chr03   9014679-9102228      * |           2
# 2    Chr03   9421216-9510839      * |           2
# 2    Chr03   9726090-9779241      * |           2
# 2    Chr03 15364219-15425632      * |           2
# 2    Chr03 15496939-15646841      * |           2
# 2    Chr03 15698107-15842550      * |           2
# 2    Chr03 15923365-16003621      * |           2
# 2    Chr03 18614795-18675274      * |           2
# -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths

####### Define the genomic ranges of the gdd bins
gddGR <- GRanges(CHR,IRanges(gdd2$POS-2.5e4,gdd2$POS+2.5e4),POS=gdd2$POS, P=gdd2$P, outlier=gdd2$outlier)

gddGRout <- unlist(reduce(split(gddGR, ~outlier)))
gddGRout$outlier <- names(gddGRout)
gddGRCand <- subset(gddGRout, outlier==2)

gddGRCand # Print the candidate regions


# GRanges object with 11 ranges and 1 metadata column:
#   seqnames            ranges strand |     outlier
# <Rle>         <IRanges>  <Rle> | <character>
#   2    Chr03     650925-840564      * |           2
# 2    Chr03     872928-979149      * |           2
# 2    Chr03   1026641-1137936      * |           2
# 2    Chr03   1173242-1369388      * |           2
# 2    Chr03   1384617-1439553      * |           2
# 2    Chr03   1537762-1596909      * |           2
# 2    Chr03   1606650-1691794      * |           2
# 2    Chr03   5507324-5567848      * |           2
# 2    Chr03   6055488-6113701      * |           2
# 2    Chr03   8796685-8846911      * |           2
# 2    Chr03 12860644-12911486      * |           2
# -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths

###### Define the genomic ranges of the ffreeze bins
ffreezeGR <- GRanges(CHR,IRanges(ffreeze2$POS-2.5e4,ffreeze2$POS+2.5e4),POS=ffreeze2$POS, P=ffreeze2$P, outlier=ffreeze2$outlier)

ffreezeGRout <- unlist(reduce(split(ffreezeGR, ~outlier)))
ffreezeGRout$outlier <- names(ffreezeGRout)
ffreezeGRCand <- subset(ffreezeGRout, outlier==2)

ffreezeGRCand # Print the candidate regions

# GRanges object with 8 ranges and 1 metadata column:
#   seqnames            ranges strand |     outlier
# <Rle>         <IRanges>  <Rle> | <character>
#   2    Chr03   5506584-5569612      * |           2
# 2    Chr03   5674399-5732422      * |           2
# 2    Chr03   6055488-6113701      * |           2
# 2    Chr03   6642576-6706657      * |           2
# 2    Chr03   7115639-7273491      * |           2
# 2    Chr03 13398490-13458285      * |           2
# 2    Chr03 15738573-15883600      * |           2
# 2    Chr03 21169441-21224799      * |           2
# -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths




#Define Overlap:  ffreeze ~ Bud Set
overlap_BF_fFReeze <- subsetByOverlaps(budflushGRCand, ffreezeGRCand)
length(overlap_BF_fFReeze)

overlap_BF_fFReeze # Print the overlapping regions

# GRanges object with 0 ranges and 1 metadata column:
#   seqnames    ranges strand |     outlier
# <Rle> <IRanges>  <Rle> | <character>
#   -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths




####Define Overlap:  gdd ~ Bud Set
overlap_BF_gdd <- subsetByOverlaps(budflushGRCand, gddGRCand)
length(overlap_BF_gdd)

overlap_BF_gdd # Print the overlapping regions

# GRanges object with 0 ranges and 1 metadata column:
#   seqnames    ranges strand |     outlier
# <Rle> <IRanges>  <Rle> | <character>
#   -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths



#####Define Overlap:  dd0 ~ Bud Set
overlap_BF_dd0 <- subsetByOverlaps(budflushGRCand, dd0GRCand)
length(overlap_BF_dd0)

overlap_BF_dd0 # Print the overlapping regions

# GRanges object with 2 ranges and 1 metadata column:
#   seqnames            ranges strand |     outlier
#       <Rle>       <IRanges>   <Rle> | <character>
# 2    Chr03   9031506-9100590      * |           2
# 2    Chr03 18515267-18776502      * |           2
# -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths



#########----------------------- Find what genes are in the overlap
# Import the GFF annotation file and make a transcript database
library(GenomicFeatures)

txdb <- makeTxDbFromGFF("Ptrichocarpa_533_v4.1.gene.gff3.gz", format="gff3") #making a transcript database

txdb

# How many chromosomes are present?
head(seqlevels(txdb))
# "Chr01" "Chr02" "Chr03" "Chr04" "Chr05" "Chr06"

# CHR="Chr03"
# Subset the database for just your chromosome of interest
seqlevels(txdb) <- CHR # subset for just your chromosome Chr03

# Reduce the transcript database to just the non-redundant gene names, instead of multiple entries for all the variant transcript types per gene
genes <- unlist(reduce(transcriptsBy(txdb, by="gene"))) 
genes$geneID <- names(genes)
#what are these genes and what do they do?? 


#### Now we'll use GenomicRanges() just like before to find the genes that overlap in the intervals of our candidate regions of overlapping lowFst/highRAiSD

BF_dd0_candGenes <- subsetByOverlaps(genes, overlap_BF_dd0)

write.table(BF_dd0_candGenes$geneID, paste0("Chr03 data/candGenes",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")

#coppy genes and go to  Popgenie.org.
head(BF_dd0_candGenes)

# GRanges object with 6 ranges and 1 metadata column:
#   seqnames            ranges strand |           geneID
# <Rle>         <IRanges>  <Rle> |      <character>
#   Potri.003G063200    Chr03   9033170-9038941      + | Potri.003G063200
# Potri.003G063300    Chr03   9042586-9051644      + | Potri.003G063300
# Potri.003G063400    Chr03   9074165-9076039      + | Potri.003G063400
# Potri.003G063500    Chr03   9100473-9105763      - | Potri.003G063500
# Potri.003G178800    Chr03 18515805-18520648      + | Potri.003G178800
# Potri.003G178900    Chr03 18521450-18528986      - | Potri.003G178900


#------------------- genes associated wtih each trait -------------

#Bud Flush
budflush_candGenes <- subsetByOverlaps(genes, budflushGRCand)

write.table(budflush_candGenes$geneID, paste0("Chr03 data/candGenes",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")

#coppy genes and go to  Popgenie.org.
head(budflush_candGenes)

# Potri.003G063200
# Potri.003G063300
# Potri.003G063400
# Potri.003G063500
# Potri.003G093700
# Potri.003G093800

#Final Freeze
ffreeze_candGenes <- subsetByOverlaps(genes, ffreezeGRCand)

write.table(ffreeze_candGenes$geneID, paste0("Chr03 data/candGenes",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")

#coppy genes and go to  Popgenie.org.
head(ffreeze_candGenes)

# Potri.003G046700
# Potri.003G047151
# Potri.003G048800
# Potri.003G048900
# Potri.003G049100
# Potri.003G049300


#dd0
dd0_candGenes <- subsetByOverlaps(genes, dd0GRCand)

write.table(dd0_candGenes$geneID, paste0("Chr03 data/candGenes",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")

#coppy genes and go to  Popgenie.org.
head(dd0_candGenes)

# Potri.003G063150
# Potri.003G063200
# Potri.003G063300
# Potri.003G063400
# Potri.003G063500
# Potri.003G066900

#gdd

gdd_candGenes <- subsetByOverlaps(genes, gddGRCand)

write.table(gdd_candGenes$geneID, paste0("Chr03 data/candGenes",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")

#coppy genes and go to  Popgenie.org.
head(gdd_candGenes)

# Potri.003G005900
# Potri.003G006000
# Potri.003G006100
# Potri.003G006200
# Potri.003G006300
# Potri.003G006500




