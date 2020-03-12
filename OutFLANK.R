library(OutFLANK)
library(vcfR)
library(data.table)
library(ggplot2)
library(qqman)
library(stringr)

setwd("~/work/Honours/")

## PREPARE INPUT DATA ##

#groups = fread("grouping_names.txt", header = FALSE)$V1
groups3 <- c("Anatolia_EF", "Steppe_wolf", "WHG")
famfile = fread("outflanktest2.fam", header=FALSE)
colnames(famfile) = c("group", "ind", paste("V", 3:6))

famfile$group_ind = paste(famfile$group, famfile$ind, sep="_")

setwd("~/work/Honours/OutFLANK_FST")

## VCF to OutFLANK format 
obj.vcfR <- read.vcfR("~/work/Honours/outflanktest2.vcf.gz")

# Character matrix containing genotypes
geno <- extract.gt(obj.vcfR) 
# Positions in bp
position <- getPOS(obj.vcfR) 
# # Chromosome information
chromosome <- getCHROM(obj.vcfR) 


save(geno, file = "OUTFLANKgenotest2.RData")
load("OUTFLANKgenotest2.RData")
 
 geno[geno %in% c("0/0", "0|0")] <- 0
 if (geno[geno  %in% c("0/1", "1/0", "1|0", "0|1")]) stop("Data contains heterozygotes!")
 geno[geno %in% c("1/1", "1|1")] <- 1

#nrowgeno pop1 8.5 hrs 100,000 three pop 35 min
countAlleles = function(famfile, geno, groups, chromosome, position) {
  Allcounts <- NULL
  start_time <- Sys.time()
  sites <-row.names(geno)
  for (i in groups){
    names <- famfile[group==i, group_ind]
    print(i)
    temp <- geno[,names]
    #for (k in 1:nrow(geno)){
    for (k in 1:100000){
      count_00 <- sum(temp[k,]=="0", na.rm=TRUE)
      count_11 <- sum(temp[k,]=="1", na.rm=TRUE)
      # count_00 <- sum(temp[k,]=="0/0", na.rm=TRUE)
      # count_11 <- sum(temp[k,]=="1/1", na.rm=TRUE)
      # count_12 <- sum(temp[k,]==c("0/1", "1/0"), na.rm=TRUE)
      # if (count_12 != 0) stop("Data contains heterozygotes!")
      alcount <- data.table(Pop = i, Chrom = chromosome[k], Position = position[k], SNP = sites[k], "0" = count_00, "1" = count_11)
      Allcounts = rbind(Allcounts, alcount)
    } 
    end_time <- Sys.time()
    runtime = end_time - start_time
    print(runtime)
  }
  return(Allcounts)
}

Allcounts = countAlleles(famfile=famfile, geno=geno, groups=groups3, chromosome = chromosome, position = position)

save(Allcounts, file = "outflanktest2_Allcount.RData")
load("outflanktest2_Allcount.RData")

 outTemp = NULL
MakeHaploidFSTMat = function (SNPmat, popNames) {  
  locusNames = unique(Allcounts$SNP)
  writeLines("Calculating FSTs, may take a few minutes...")
  for (i in 1:length(locusNames)) {
    SNPSITE<- SNPmat[SNP==locusNames[i]]
    SNPSITE[, SNP := NULL]
    SNPSITE[, Pop := NULL]
    SNPSITE[, Chrom := NULL]
    SNPSITE[, Position:= NULL]
    input <- as.matrix(SNPSITE)
    if (dim(input)[1] != length(popNames)) {
    print("Error: your population names do not match your SNP matrix")
    break
    }
      Nocorr <- unlist(WC_FST_FiniteSample_Haploids_2AllelesB_NoSamplingCorrection(input))
      Corr <- unlist(WC_FST_FiniteSample_Haploids_2AllelesB_MCW(input))
      out <- data.table(Pop1 = popNames[1], Pop2 = popNames[2], Pop3 = popNames[3], 
                        LocusName = locusNames[i], He = Corr[1], FST = Corr[3], 
                        T1 = Corr[4], T2 = Corr[5], FSTNoCorr = Nocorr[3], 
                        T1NoCorr = Nocorr[4], T2NoCorr = Nocorr[5], meanAlleleFreq = Corr[2], 
                        Chromosome = SNPmat$Chrom[i], Position = SNPmat$Position[i])
      outTemp = rbind(outTemp, out)
    if (i%%10000 == 0) {
      print(paste(i, "done of", length(locusNames)))
      }
  }
  save(outTemp, file= "outTemp.Rdata")
  return(outTemp)
}

#load("outTemp.Rdata")

my_fst <- MakeHaploidFSTMat(SNPmat = Allcounts, popNames = groups3)

# OutFLANKHaploid = function(SNPmat, FSTOutFLANK){}
Allcounts[, sampleSize := `0` + `1`]
sampleSizes = dcast(Allcounts, SNP ~ Pop, value.var = "sampleSize")
colnames(sampleSizes) = c("LocusName", "pop1", "pop2", "pop3")

my_fst = merge(my_fst, sampleSizes, by = "LocusName", all.x=TRUE)
my_fst[,totalSampleSize := pop1 + pop2 + pop3]
my_fst[,minSampleSize := pmin(pop1, pop2, pop3, na.rm=T)]

## Data check plots ##
ggplot(my_fst, aes(minSampleSize, FSTNoCorr)) + geom_point() + geom_smooth()

ggplot(my_fst[minSampleSize >= 2 & He > 0.1], aes(FST, FSTNoCorr, col=minSampleSize)) + 
  geom_point() + geom_smooth() + geom_abline(aes(intercept=0, slope=1))

ggplot(my_fst[minSampleSize >= 2], aes(He, FST, col=minSampleSize)) + 
  geom_point() + geom_smooth() 

ggplot(my_fst[totalSampleSize > 20], aes(FST, FSTNoCorr, col=totalSampleSize)) + 
  geom_point() + geom_smooth() + geom_abline(aes(intercept=0, slope=1))

ggplot(my_fst[ID_filter], aes(FST, FSTNoCorr, col=totalSampleSize)) + 
  geom_point() + geom_smooth() + geom_abline(aes(intercept=0, slope=1))

## FILTER ##
ID_filter = my_fst$minSampleSize >= 2 & !is.nan(my_fst$FST) & !is.nan(my_fst$FSTNoCorr) & my_fst$FSTNoCorr > 0 

out_trim <- OutFLANK(my_fst[ID_filter,.(LocusName, FST, T1, T2, FSTNoCorr, T1NoCorr, T2NoCorr, He)], NumberOfSamples=3, qthreshold = 0.05, Hmin = 0.1, LeftTrimFraction = 0.05, RightTrimFraction = 0.05)
str(out_trim)

## CHECK FIT ##
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
#right tail
OutFLANKResultsPlotter(out_trim , withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         TRUE, RightZoomFraction = 0.15, titletext = NULL)
#p-value hist
hist(out_trim$results$pvaluesRightTail, n=100)

#p-val for all loci
P1 <- pOutlierFinderChiSqNoCorr(my_fst[ID_filter,], Fstbar = out_trim$FSTNoCorrbar, 
                                dfInferred = out_trim$dfInferred, qthreshold = 0.1, Hmin=0.1)
head(P1)

#OUTPUT#
 #check output
my_out <- P1$OutlierFlag==TRUE
plot(P1$He, P1$FST, pch=19, col=rgb(0,0,0,0.1))
points(P1$He[my_out], P1$FST[my_out], col="blue")

hist(P1$pvaluesRightTail)

#highlight outliers
plot(P1$LocusName[P1$He>0.1], P1$FST[P1$He>0.1],
     xlab="Position", ylab="FST", col=rgb(0,0,0,0.2))
points(P1$LocusName[my_out], P1$FST[my_out], col="magenta", pch=20)


manhattan(P1[my_out],chr="Chromosome",bp="Position", snp = "LocusName",p="pvalues",
                            col="black", logp=TRUE,ylab="p-value",
                            highlight = NULL)

# manhattan(P1[my_out],chr="Chromosome",bp="Position", snp = "LocusName",p="FSTNoCorr",
#                             col="black", logp=FALSE,ylab="FST",
#                             highlight = NULL)

## Highlight OAGR SNPS ##
# OAGR_SNPS <- read_excel("~/work/Honours/FST/OAGR_SNPS.xlsx")
# OAGR_SNPS$END <- OAGR_SNPS$END*1e3
# OAGR_SNPS$START <- OAGR_SNPS$START*1e3
# 
# OAGRSNP_LIST <- NULL
# for (i in 1:length(P1[my_out])){
#   for (j in 1:length(OAGR_SNPS$CHROM)) {
#     if (P1$Position[i] >= OAGR_SNPS$START[j] & P1$Position[i] <= OAGR_SNPS$END[j]){
#       if (P1$Chromosome[i] == OAGR_SNPS$CHROM[j]) {
#         OAGRSNP_LIST <- list.append(OAGRSNP_LIST, toString(P1$LocusName[i]))
#       }
#     }
#   }
# }
# 
# for (i in 1:length(OAGRSNP_LIST)){
#   if (OAGRSNP_LIST[i] %in% P1$LocusName){
#   } else{
#     print("Missing values")
#   }
# }
# manhattan_plot <- manhattan(P1[my_out],chr="Chromosome",bp="Position", snp = "LocusName",p="pvalues",
#                             col="black", logp=TRUE,ylab="p-value",
#                             highlight = OAGRSNP_LIST)