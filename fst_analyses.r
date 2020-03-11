
library(data.table)
library(ggplot2)
library(stringr)
library(stringi)
library(foreach)
library(matrixStats)
library(qvalue)
library(IRanges)


##--------------------------------------------------##
#Haploid Fst estimator
##--------------------------------------------------##

fst_hap <- function(dat, pops){
   dt.pop <- dat[Population %in% pops]
   r <- length(pops) #number of populations

   #make data matrices
   d.n <- dcast(dt.pop, SNP~Population, value.var="Count")
   d.N <- dcast(dt.pop, SNP~Population, value.var="N")
   d.f <- dcast(dt.pop, SNP~Population, value.var="Freq")

   #perform calculations
   n_sum <- rowSums(d.N[, 2:(r+1), with=F])
   n_ave <- n_sum/r
   n_sum2 <- rowSums(d.N[, (.SD)^2, .SDcols=2:(r+1)])
   n_c <- (n_sum - n_sum2/n_sum)/(r-1)

   p_ave <- rowSums(d.n[, 2:(r+1), with=F])/n_sum

   s2 <- rowSums(d.N[, .SD, .SDcols=2:(r+1)] *
                 d.f[, ((.SD)-p_ave)^2, .SDcols=2:(r+1)])/((r-1)*n_ave)

   He <- 2*p_ave*(1-p_ave)

   T1.corr <- s2 - 1/(n_ave-1)*(p_ave*(1-p_ave)-(s2*(r-1)/r))
   T2.corr <- (n_c-1)*(p_ave*(1-p_ave))/(n_ave-1) + (1+(r-1)*(n_ave-n_c)/(n_ave-1))*s2/r
   FST.corr <- T1.corr/T2.corr

   T1 <- s2
   T2 <- s2/r + (p_ave*(1-p_ave))
   FST <- T1/T2

   setnames(d.n, c("X",  paste0("COUNT", 1:length(pops))))
   setnames(d.N, c("X",  paste0("TOTAL", 1:length(pops))))
   setnames(d.f, c("X",  paste0("f", 1:length(pops))))

   dt.out <- foreach(p=pops) %do% rep(p, dt.pop[, length(unique(SNP))])
   setDT(dt.out)
   setnames(dt.out, paste0("Pop", 1:length(pops)))

   data.table(dt.out, dt.pop[Population==pops[1], .(Chr, Pos, SNP)],
              d.n[, 2:(r+1), with=F], d.N[, 2:(r+1), with=F],
              d.f[, 2:(r+1), with=F], He, T1, T2, FST, T1.corr,
              T2.corr, FST.corr)[order(Chr, Pos)]
}

##--------------------------------------------------##
#PBS estimator
##--------------------------------------------------##

#calculate pairwise fst for all combinations of listed populations

#the populations is a list of vectors with the names all trios used in estimation
#focal group must be listed first
pbs.est <- function(dat, pops, window=200000, step=100000, min.N=5, min.He=0.01){

   #loop over all lists of trios
   completed.duos <- c() #keep record of previously computed duos

   #calculate all unique starting points for windows
   unq.start <- seq(0, window-1, step)

   #create window labels
   window.labs <- foreach(US=unq.start) %do% {
      i <- which(unq.start==US)
      windows <- unique(dat[, c(0, seq(US, max(Pos)+window, window))])
      start.labs <- windows[1:(length(windows)-1)]
      end.labs <- windows[2:length(windows)]-1
      list(WINDOW=windows,
           POS=data.table(I=paste(i, 1:length(start.labs), sep="."),
                          Start=start.labs,
                          End=end.labs))
   }
   names(window.labs) <- unq.start

   #keep track of all FST calculations
   T.mat <- NULL

   #estimate pbs
   pbs.dt <- foreach(k=pops, .combine=rbind) %do% {

      pbs.trio=paste(k, collapse=",")
      all.fst.duos <- combn(k, 2, simplify=F)

      print(paste("calculating PBS for", pbs.trio))

      trio.now.mat <- NULL
      cc=0
      #estimate FST in bins, loop over different starting points
      for(DUO in all.fst.duos){

         duo.sort <- paste0(sort(DUO), collapse=",")
         print(paste("calculating FST for duo", duo.sort))

         cc=cc+1
         if(duo.sort %in% completed.duos){
            t.mat <- T.mat[Fst.Duo==duo.sort]
            print(paste0("Reusing ", duo.sort))
         }else{
            fst.now <- fst_hap(dat=dat, pops=DUO)[TOTAL1>min.N & TOTAL2>min.N]
            #loop over different window steps
            t.mat <- foreach(US=unq.start, .combine=rbind) %do% {
               print(paste("step start:", US))
               windows.now <- window.labs[[which(names(window.labs)==US)]]
               fst.now[, WINDOW:=cut(Pos, windows.now$WINDOW,
                                     labels=windows.now$POS$I), by="Chr"]
               out <- fst.now[, .(Fst=sum(T1)/sum(T2),
                                  n.SNP=.N), by=c("Chr", "WINDOW")]
               out[, Txy:=-log(1-Fst)]
               data.table(Fst.Duo=duo.sort, I=cc, out)
            }
            T.mat <- rbind(T.mat, t.mat)
            completed.duos <- c(completed.duos, duo.sort)
         }
         trio.now.mat <- rbind(trio.now.mat, t.mat)
      }

      EST <- function(M){
         M[, Pop.Obs:=.N, by=c("Chr", "WINDOW")]

         NN <- dcast(M[Pop.Obs==3], Chr+WINDOW~Fst.Duo, value.var="n.SNP")
         setnames(NN, c("Chr", "WINDOW", "N12", "N13", "N23"))

         MM <- dcast(M[Pop.Obs==3], Chr+WINDOW~Fst.Duo, value.var="Txy")
         setnames(MM, c("Chr", "WINDOW", "T12", "T13", "T23"))
         MM[, PBS:=(T12+T13-T23)/2]
         out <- merge(NN, MM, by=c("Chr", "WINDOW"))
         out
      }

      XX <- EST(trio.now.mat)
      #add start and end for each region
      YY <- foreach(i=window.labs, .combine=rbind) %do% i$POS

      out <- merge(XX, YY, by.x="WINDOW", by.y="I")

      data.table(Trio=pbs.trio, out)[, .(Trio, Chr, Start, End, WINDOW,
                                         N12, N13, N23, T12, T13, T23,
                                         PBS)][order(Chr, Start)]
   }
   print("finished calculating PBS")
   return(pbs.dt)
}


##--------------------------------------------------##
#outflank functions
##--------------------------------------------------##

FstDistPlotter = function(df, FSTlist, FSTbar, binwidth=0.005, titletext=NULL){
   xPlotUpperBound=ceiling(max(FSTlist)*100)/100
   breakslist=seq(0,xPlotUpperBound+binwidth,by=binwidth)
   breaks = length(breakslist)

   x = breakslist
   y=rep(0,length(x))
   for(i in 1:breaks) y[i] = pchisq(((i-.5)*binwidth)/FSTbar*df , df=df) - pchisq((((i-1.5)*binwidth))/FSTbar*df , df=df)
   y=length(FSTlist)*y

   hist(FSTlist,col="darkgoldenrod1", breaks=breakslist, prob=F, xlab="Fst",  main=titletext)

   lines(x,y,col="darkblue", lwd=3)
}

EffectiveNumberSamplesMLE=function(FST, FSTbar, Npops){

   localNLLAllData <- function(dfInferred){
      localNLLOneLocus <- function(FST){
         negLLdfFstTrim(FST, dfInferred, FSTbar)
      }
      sum(localNLLOneLocus(FST), na.rm=T)
   }

   #optim(Npops, localNLLAllData, lower=0, method="L-BFGS-B")$par
   optim(Npops, localNLLAllData, method="BFGS")$par
}

IncompleteGammaFunction=function(a, z) {
   #equivalence to Mathematica Gamma[a,z] according to
   #   http://r.789695.n4.nabble.com/Incomplete-Gamma-function-td833545.html
   pgamma(z, a, lower=FALSE)*gamma(a)
}

negLLdfFstTrim=function(FST, dfInferred, FSTbar){
   df <- dfInferred
   HighTrimPoint <- max(FST)
   LowTrimPoint <- min(FST)

   1/(2*FSTbar) * (df*FST+df*FSTbar*log(2)-df*FSTbar*log(df) -
                      (df-2)*FSTbar*log(FST) +
                      df*FSTbar*log(FSTbar) +
                      2*FSTbar*log(-IncompleteGammaFunction(df/2, df*HighTrimPoint/(2*FSTbar)) +
                                      IncompleteGammaFunction(df/2, df*LowTrimPoint/(2*FSTbar))))
}

pTwoSidedFromChiSq=function(x, df){
   #pOneSided=pchisq(x, df)
   #ifelse(pOneSided>.5, (1-pOneSided)*2, pOneSided*2)
   pOneSided=1-pchisq(x, df)
   ifelse(pOneSided>.5, NA, pOneSided*2)
}

pOutlierFinderChiSq <- function(dt.now, FSTbar, dfInferred){
   dt.now[KEEP==T, p2sided:=pTwoSidedFromChiSq(FST*(dfInferred)/FSTbar, dfInferred)]
   dt.now[KEEP==T, pRightTail:=1-pchisq(FST*(dfInferred)/FSTbar, dfInferred)]
   #dt.now[KEEP==T, qRightTail.BH:=p.adjust(pRightTail, method="BH")]
   dt.now[KEEP==T, q2sided:=qvalue(p2sided, pi0.method="bootstrap")$qvalues]
   dt.now[KEEP==T, qRightTail:=qvalue(pRightTail, pi0.method="bootstrap")$qvalues]
}

OutFLANK <- function(dd, pops.now, Hmin=0.1, Nthreshold=5, qthreshold=0.05,
                     LeftTrimFraction=0.025, RightTrimFraction=0.025){

   #subset datatable
   pop.cols <- colnames(dd)[grep("^Pop", colnames(dd))]
   ff <- foreach(pp=pop.cols, .combine='+') %do% dd[, get(pp) %in% pops.now]
   dt.now <- dd[which(ff==length(pops.now))]

   #remove NAs, SNPs with low counts and low He
   dt.now[, OutlierFlag:=ifelse(is.na(FST), NA, F)]
   dt.now[, KEEP:=!(is.na(FST) | He<Hmin)]

   N.cols <- colnames(dt.now)[grep("^TOTAL", colnames(dt.now))]
   nn <- foreach(NN=N.cols, .combine='+') %do% dt.now[, get(NN)<Nthreshold]
   dt.now[which(nn>0), KEEP:=FALSE]

   NLoci <- dt.now[KEEP==T, .N]

   #determine set of trimmed SNPs
   TrimPoints <- dt.now[KEEP==T, quantile(FST, c(LeftTrimFraction,
                                                 1-RightTrimFraction))]
   dt.now[, TRIM.KEEP:=KEEP]
   dt.now[FST<TrimPoints[1] | FST>TrimPoints[2], TRIM.KEEP:=F]

   #finding dfInferred and Fstbar iteratively
   dt.now[, OutlierFlag:=F]

   keepGoing=TRUE
   count = 0

   #keep track of changing parameters
   #trimpoints <- list()
   trimpoints <- TrimPoints
   df <- c()
   fstbar <- c()

   while(keepGoing){
      count=count+1
      print(paste("Iteration", count))
      if(count>19) {
         keepGoing=FALSE
         writeLines("Exceeded iteration maximum.") ###Try with increased maximum value for count two lines above.
      }

      FSTbar <- dt.now[TRIM.KEEP==T, mean(FST, na.rm=T)]
      dfInferred <- EffectiveNumberSamplesMLE(dt.now[TRIM.KEEP==T, FST],
                                              FSTbar, Npops=length(pops.now))

      pOutlierFinderChiSq(dt.now, FSTbar, dfInferred)

      #determine if any new outliers
      any.new.ol <- dt.now[OutlierFlag==F & qRightTail<qthreshold, .N]

      if(any.new.ol>0){
         #update parameters
         #dt.now[qRightTail<qthreshold, OutlierFlag:=T]
         dt.now[, OutlierFlag:=ifelse(qRightTail<qthreshold, T, F)]
         dt.now[, KEEP:=!(is.na(FST) | He<Hmin | OutlierFlag==T)]

         NLoci <- dt.now[KEEP==T, .N]

         TrimPoints <- dt.now[KEEP==T, quantile(FST, c(LeftTrimFraction,
                                                       1-RightTrimFraction))]
         dt.now[, TRIM.KEEP:=KEEP]
         dt.now[FST<TrimPoints[1] | FST>TrimPoints[2], TRIM.KEEP:=F]

         trimpoints <- rbind(trimpoints, TrimPoints)
         df <- c(df, dfInferred)
         fstbar <- c(fstbar, FSTbar)
      }else{
         keepGoing=FALSE
      }
   }

   if(count>19) writeLines("Loop iteration limit exceeded.")

   list(DT=dt.now, DF=df, FSTbar=fstbar, TrimPoints=trimpoints)
}

#source("~/Dropbox/oagr_manuscript/WorkDir/R/outflank.r")

##--------------------------------------------------##
#Read in fst data
##--------------------------------------------------##

setwd("~/Dropbox/oagr_manuscript/WorkDir/genotypeData/FINAL_Data/Freq_files")
#system("plink2 --bfile ../ancient_west_eurasians --freq --family --remove-fam fam_remove_subset.txt")

#read in allele frequency data
bed <- fread("plink.frq.strat")

#create input data.table and set YRI as outgroup
dt <- bed[, .(Chr=CHR, SNP, Population=CLST, Freq=MAF,
              Count=as.double(MAC), N=as.double(NCHROBS))]

#reset counts and frequencies for pseudohaploid populations (IGNORE warnings)
dt[!(Population %in% c("CHB", "CEU", "TSI", "FIN", "YRI")), Count:=Count/2]
dt[!(Population %in% c("CHB", "CEU", "TSI", "FIN", "YRI")), N:=N/2]
dt[!(Population %in% c("CHB", "CEU", "TSI", "FIN", "YRI")), Freq:=Count/N]

get.pos <- data.table(SNP=unique(dt$SNP),
                      Pos=as.numeric(sapply(stri_split(unique(dt$SNP), fixed="_"), "[", 2)))

dt <- merge(dt, get.pos, by="SNP")

rm(bed)

##--------------------------------------------------##
# outflank estimation
##--------------------------------------------------##

sfs <- fread("~/Dropbox/oagr_manuscript/WorkDir/PolySel_Outputs/Plots/SI/effective.N.info_FINAL.txt")
sfs[Population=="Levant", Population:="Levant_EF"]
small <- sfs[eff.N<10, Population]

#run two estimation steps, with and without YRI
dt.list <- foreach(moderns=list(c("CEU", "FIN", "TSI", "CHB", "YRI"),
                    c("CEU", "FIN", "TSI", "CHB"))) %do% {

   all.pops <- dt[!(Population %in% c(moderns, small)), unique(Population)]

   dt.fst <- fst_hap(dt, pops=all.pops)

   #run outflank
   print(paste("running outflank for", paste(all.pops, collapse=", ")))
   OutFLANK(dd=dt.fst, pops.now=all.pops, qthreshold=0.05,
            Nthreshold=5, Hmin=0.1,
            LeftTrimFraction=0.025, RightTrimFraction=0.025)
}
names(dt.list) <- c("YRI excluded", "YRI included")

save(dt.list, file="~/Dropbox/oagr_manuscript/WorkDir/FST.results.rdata")

##--------------------------------------------------##
# compare with sweep windows
##--------------------------------------------------##

#load("~/Dropbox/oagr_manuscript/WorkDir/FST.results.rdata")

sweeps <- fread("/Users/raytobler/Dropbox/oagr_manuscript/WorkDir/SWEEP.LIST.TOP17.pathways.FINAL.txt")
ss <- strsplit(sweeps$IGV.window, ":")
ss2 <- strsplit(sapply(ss, "[", 2), "-")
sweeps[, Chr:=as.numeric(sapply(ss, "[", 1))]
sweeps[, Start:=as.numeric(sapply(ss2, "[", 1))]
sweeps[, End:=as.numeric(sapply(ss2, "[", 2))]

sw <- sweeps[, .(Chr=unique(Chr), Start=unique(Start), End=unique(End)), by="SweepID"]

#get all SNPs overlapping sweeps
gg.dt <- foreach(k=1:length(dt.list), .combine=rbind) %do% {
   DT <- data.table(Type=names(dt.list)[k],
                    dt.list[[k]]$DT[, .(SNP, He, FST,
                                        OutlierFlag, KEEP,
                                        pRightTail, qRightTail)])

   ss <- stri_split(DT$SNP, fixed="_")
   DT[, Chr:=sapply(ss, "[", 1)]
   DT[, Pos:=as.numeric(sapply(ss, "[", 2))]

   sw.f <- foreach(chr.now=1:22, .combine=rbind) %do% {
      print(paste("Generating for Chr", chr.now))
      sw.now <- sw[Chr==chr.now]
      DT.now <- DT[Chr==chr.now]

      subject <- IRanges(DT.now$Pos, width=rep(1, nrow(DT.now)))
      query <- IRanges(sw.now$Start, width=sw.now[, End-Start])

      ol <- findOverlaps(query, subject)

      data.table(DT.now[ol@to], sw.now[ol@from, .(SweepID)])
   }
   out <- rbind(sw.f,
                data.table(DT[!(SNP %in% sw.f$SNP)],
                           SweepID="Background"))

   out
}

#SNPs not in sweeps
mm <- gg.dt[, .(MM=mean(-log10(pRightTail), na.rm=T)),
            by=c("Type", "SweepID")]

new.ord <- mm[Type=="YRI included"][order(MM), SweepID]
gg.dt[, SweepID:=factor(SweepID, new.ord)]
gg.dt[, COL:=factor(ifelse(SweepID=="Background", 1, 0))]

gg.dt[, neglog10p:=-log10(pRightTail)]
quid <- gg.dt[!is.na(neglog10p), .(Lower=quantile(neglog10p, 0.25),
                                   Upper=quantile(neglog10p, 0.75),
                                   Median=median(neglog10p)), by=c("SweepID", "Type")]

#rank test on p-values vs background
rank.test <- foreach(sw.now=sw$SweepID, .combine=rbind) %do% {
   print(sw.now)
   yi.sw <- gg.dt[Type=="YRI included" & SweepID==sw.now & !(is.na(pRightTail))]
   yx.sw <- gg.dt[Type=="YRI excluded" & SweepID==sw.now & !(is.na(pRightTail))]
   yi.bk <- gg.dt[Type=="YRI included" & SweepID=="Background"  & !(is.na(pRightTail))]
   yx.bk <- gg.dt[Type=="YRI excluded" & SweepID=="Background"  & !(is.na(pRightTail))]

   wt1 <- wilcox.test(yi.sw[, neglog10p], yi.bk[, neglog10p],
                      paired=F, alternative="greater")
   wt2 <- wilcox.test(yx.sw[, neglog10p], yx.bk[, neglog10p],
                      paired=F, alternative="greater")
   wt3 <- wilcox.test(yi.sw[, neglog10p], yx.sw[, neglog10p],
                      paired=F, alternative="greater")

   wt4 <- wilcox.test(yi.sw[, FST], yi.bk[, FST],
                      paired=F, alternative="greater")
   wt5 <- wilcox.test(yx.sw[, FST], yx.bk[, FST],
                      paired=F, alternative="greater")
   wt6 <- wilcox.test(yi.sw[, FST], yx.sw[, FST],
                      paired=F, alternative="greater")

   data.table(SweepID=sw.now, P_YRI=wt1$p.value, P_noYRI=wt2$p.value, P_vs=wt3$p.value,
              Fst_YRI=wt4$p.value, Fst_noYRI=wt5$p.value, Fst_vs=wt6$p.value)
}




rank.test[, SvsB:="Neither"]
rank.test[P_YRI<0.05 & P_noYRI<0.05, SvsB:="Both"]
rank.test[P_YRI<0.05 & P_noYRI>=0.05, SvsB:="YRI present only"]
rank.test[P_YRI>=0.05 & P_noYRI<0.05, SvsB:="YRI absent only"]

rank.test[, YPvsYA:="Non-sig"]
rank.test[P_vs<0.05, YPvsYA:="Sig"]

squid <- merge(merge(quid[Type=="YRI excluded", .(SweepID, Median, Lower, Upper)],
                     quid[Type=="YRI included", .(SweepID, Median, Lower, Upper)],
                     by="SweepID"), rank.test, by="SweepID", all=T)

squid[, `Sweep vs Backgound`:=paste0(SvsB, " [", .N, "]"), by="SvsB"]
squid[, `YRI present vs YRI absent`:=paste0(YPvsYA, " [", .N, "]"), by="YPvsYA"]

squid[is.na(SvsB), `Sweep vs Backgound`:="Background"]
squid[is.na(YPvsYA), `YRI present vs YRI absent`:="Background"]

squid[, `Sweep vs Backgound`:=factor(`Sweep vs Backgound`,
                                     levels=squid[, sort(unique(`Sweep vs Backgound`))[c(1,3,2,5,4)]])]

squad <-  merge(squid,
                fread("~/Dropbox/oagr_manuscript/WorkDir/sweep.classification.txt")[, .(SweepID, Earliest_evidence)],
                by="SweepID")


ggplot(data=squad,
       aes(x=Median.x, y=Median.y, shape=`Sweep vs Backgound`,
           col=Earliest_evidence)) +
   geom_abline(intercept=0, slope=1, size=0.75) +
   geom_hline(yintercept=squid[SweepID=="Background", Median.y],
              color="black", size=0.25, linetype=2) +
   geom_vline(xintercept=squid[SweepID=="Background", Median.x],
              color="darkred", size=0.25, linetype=2) +
   #geom_segment(aes(x=Lower.x, y=Median.y, xend=Upper.x ,
   #                 yend=Median.y, colour=`Sweep vs Backgound`),
   #             size=0.15) +
   #geom_segment(aes(x=Median.x, y=Lower.y, xend=Median.x ,
   #                 yend=Upper.y, colour=`Sweep vs Backgound`),
   #             size=0.15) +
   scale_colour_manual(values=c("black", "gray50", "red",
                                "orange", "blue", "green")) +
   geom_point(alpha=0.5, size=2) +
   theme_bw() +
   xlab("-log10(p) without YRI") + ylab("-log10(p) with YRI") +
   theme(legend.position="bottom",
         legend.box="vertical",
         legend.spacing=unit(0.05, 'cm'))

squib <- squad[, .N, by=.(SvsB, Earliest_evidence)]

squib[, Prop:=N/sum(N), by="Earliest_evidence"]

ggplot(squib,
       aes(y=Prop, x=Earliest_evidence, fill=SvsB)) +
   geom_bar(stat="identity") +
   scale_x_discrete(drop=FALSE)

plot.path <- "~/Dropbox/oagr_manuscript/WorkDir/PolySel_Outputs/Plots/SI"
ggsave(filename=file.path(plot.path, "FST_OutFLANK.FINAL.YRI_vs_noYRI.NEW.pdf"),
       width=8, height=8.25, dpi=300)


#just YRI
gg.dtx <- merge(gg.dt[Type=="YRI included" & KEEP==T],
                rank.test[, .(SweepID, P_YRI)],
                by="SweepID", all.x=T)

gg.dtx[SweepID=="Background", `Rank test`:="BG"]
gg.dtx[P_YRI>0.05, `Rank test`:="P >= 0.05"]
gg.dtx[P_YRI<0.05, `Rank test`:="P < 0.05"]

sw.ord <- gg.dtx[, mean(-log10(pRightTail), na.rm=T),
                 by="SweepID"][order(V1), SweepID]

gg.dtx[, SweepID.ord:=factor(SweepID, levels=sw.ord)]


ggplot(gg.dtx[Type=="YRI included" & KEEP==T],
       aes(x=SweepID.ord, y=-log10(pRightTail))) +
   #geom_boxplot() +
   geom_violin(aes(fill=`Rank test`, col=`Rank test`), alpha=0.5) +
   geom_point(data=gg.dt[OutlierFlag==T],
              aes(x=SweepID.ord, y=-log10(pRightTail)),
              col="red", alpha=0.5, size=0.5) +
   geom_hline(yintercept=gg.dtx[, mean(-log10(pRightTail), na.rm=T)]) +
   scale_fill_manual(values=c("blue", "red", "gray50")) +
   scale_color_manual(values=c("blue", "red", "gray50")) +
   stat_summary(fun.y=median, colour="gray10", geom="point",
                shape=18, size=2, show.legend=FALSE) +
   stat_summary(fun.y=mean, colour="darkred", geom="point",
                shape=18, size=2, show.legend=FALSE) +
   stat_summary(fun.data=mean_se, geom="errorbar", colour="darkred") +
   theme_bw() +
   theme(axis.text.x=element_text(angle=45, hjust=1),
         axis.title.x=element_blank(),
         legend.position="bottom",
         legend.margin=unit(0, "line"))

plot.path <- "~/Dropbox/oagr_manuscript/WorkDir/PolySel_Outputs/Plots/SI"
ggsave(filename=file.path(plot.path, "FST_OutFLANK.FINAL.NEW.pdf"),
       width=10, height=6, dpi=300)





rank.test <- foreach(sw.now=sw$SweepID, .combine=rbind) %do% {
   print(sw.now)
   yi.sw <- gg.dt[Type=="YRI included" & SweepID==sw.now & !(is.na(FST))]
   #yx.sw <- gg.dt[Type=="YRI excluded" & SweepID==sw.now & !(is.na(pRightTail))]
   yi.bk <- gg.dt[Type=="YRI included" & SweepID=="Background"  & !(is.na(FST))]
   #yx.bk <- gg.dt[Type=="YRI excluded" & SweepID=="Background"  & !(is.na(pRightTail))]

   wt1 <- wilcox.test(yi.sw[, FST], yi.bk[, FST],
                      paired=F, alternative="greater")

   data.table(SweepID=sw.now, P_YRI=wt1$p.value)
}

gg.dtx <- merge(gg.dt[Type=="YRI included"],
                rank.test[, .(SweepID, Fst_YRI)],
                by="SweepID", all.x=T)

gg.dtx[SweepID=="Background", `Rank test`:="BG"]
gg.dtx[Fst_YRI>0.05, `Rank test`:="p>=0.05"]
gg.dtx[Fst_YRI<0.05, `Rank test`:="p<0.05"]

gg.dtx[, `Fst q`:=cut(qRightTail, c(0, 0.05, 0.20, 1),
                      c("q<0.05", "0.05<=q<0.20"))]

sw.ord <- gg.dtx[Type=="YRI included", mean(FST, na.rm=T),
                 by="SweepID"][order(V1), SweepID]

gg.dtx[, SweepID.ord:=factor(SweepID, levels=sw.ord)]
#gg.dt[, SweepID.ord:=factor(SweepID, levels=sw.ord)]


ggplot(gg.dtx[Type=="YRI included"],
       aes(x=SweepID.ord, y=FST)) +
   #geom_boxplot() +
   geom_violin(aes(fill=`Rank test`, col=`Rank test`),
               alpha=0.5, width=1.1, lwd=0.25) +
   geom_boxplot(outlier.shape=NA, width=0.15, lwd=0.25) +
   geom_point(data=gg.dtx[qRightTail<0.2],
              aes(x=SweepID.ord, y=FST, size=`Fst q`),
              alpha=0.5) +
   geom_hline(yintercept=gg.dtx[, mean(FST, na.rm=T)],
              col="darkred") +
   scale_fill_manual(values=c("blue", "red", "gray75")) +
   scale_colour_manual(values=c("blue", "red", "gray75")) +
   scale_size_manual(values=c(1.25, 0.5)) +
   #stat_summary(fun.y=median, colour="gray10", geom="point",
   #             shape=18, size=3, show.legend=FALSE) +
   stat_summary(fun.y=mean, colour="darkred", geom="point",
                shape=18, size=3, show.legend=FALSE) +
   stat_summary(fun.data=mean_se, geom="errorbar",
                colour="darkred") +
   theme_bw() +
   theme(axis.text.x=element_text(angle=45, hjust=1),
         axis.title.x=element_blank(),
         legend.position="bottom")

plot.path <- "~/Dropbox/oagr_manuscript/WorkDir/PolySel_Outputs/Plots/SI"
ggsave(filename=file.path(plot.path, "FST_OutFLANK.FINAL.NEW.pdf"),
       width=10, height=6, dpi=300)

##--------------------------------------------------##
# are outliers driven by YRI
##--------------------------------------------------##

cuts <- c(0, 5, 100)/100
log.cuts <- -log10(labs)
labs <- c("p<0.05", "p>=0.5")

hibs <- merge(merge(gg.dt,
                    rbind(rank.test[, .(Type="YRI included", SweepID, p=Fst_YRI,
                                        'Wilcox Test'=cut(Fst_YRI, cuts, labs),
                                        Fst_vs, pVs=cut(Fst_vs, cuts, labs))],
                          rank.test[, .(Type="YRI excluded", SweepID, p=Fst_noYRI,
                                        'Wilcox Test'=cut(Fst_noYRI, cuts, labs),
                                        Fst_vs, pVs=cut(Fst_vs, cuts, labs))]),
                    by=c("SweepID", "Type"), all.x=T),
              fread("~/Dropbox/oagr_manuscript/WorkDir/sweep.classification.txt")[, .(SweepID, Earliest_evidence)],
              all=T, by="SweepID")

hibs[is.na(Earliest_evidence), Earliest_evidence:="BG"]
#hibs[, Upper:=quantile(FST, 0.95, na.rm=T)]
#hibs[, Lower:=quantile(FST, 0.05, na.rm=T)]
#hibs[, FST.100:=FST*100]
zz <- hibs[Earliest_evidence=="BG",
           .(Z=mean(FST, na.rm=T)), by="Type"]

sw.ord <- hibs[Type=="YRI included",
               .(MM=mean(FST, na.rm=T)),
               by="SweepID"][order(MM), SweepID]

hibs[, SweepID.ord:=factor(SweepID, levels=sw.ord)]
hibs[Earliest_evidence=="BG", `Wilcox Test`:="Background"]


ggplot(hibs,
       aes(x=SweepID.ord, y=FST, fill=`Wilcox Test`, col=`Wilcox Test`)) +
   geom_violin(alpha=0.5) +
   stat_summary(fun.y=mean, geom="point", colour="darkred",
                shape=18, size=3, show.legend=FALSE) +
   stat_summary(fun.data=mean_se, geom="errorbar", colour="darkred") +
   #geom_errorbar(data=hibs, aes(ymin=Lower, ymax=Upper)) +
   facet_grid(Type~Earliest_evidence, space="free_x", scales="free_x") +
   geom_hline(data=zz, aes(yintercept=Z), col="darkred") +
   #stat_summary(fun.data=mean_se, geom="errorbar", colour="darkred") +
   theme_bw() +
   #scale_y_log10() +
   scale_fill_manual(values=c("red", "gray50", "blue")) +
   scale_color_manual(values=c("red", "gray50", "blue")) +
   theme(axis.text.x=element_text(angle=45, hjust=1),
         axis.title.x=element_blank(),
         legend.position="bottom",
         strip.text.x=element_text(angle=90))

plot.path <- "~/Dropbox/oagr_manuscript/WorkDir/PolySel_Outputs/Plots/SI"
ggsave(filename=file.path(plot.path, "FST_OutFLANK.YRI_vs_nonYRI.pdf"),
       width=10, height=6, dpi=300)



##--------------------------------------------------##
# pbs estimation
##--------------------------------------------------##

pop.list <- list(c("WHG", "CHB", "YRI"),
                 c("Anatolia_EF", "CHB", "YRI"),
                 c("Steppe", "CHB", "YRI"),
                 c("CentralEurope_EF", "CHB", "YRI"),
                 c("CentralEurope_LNBA", "CHB", "YRI"),
                 c("CEU", "CHB", "YRI"),
                 c("Anatolia_EF", "WHG", "CHB"),
                 c("CentralEurope_EF", "WHG", "CHB"),
                 c("CentralEurope_LNBA", "WHG", "CHB"),
                 c("CEU", "WHG", "CHB"))

pbs.out <- pbs.est(dat=dt, pops=pop.list)

#rank test
rank.test <- foreach(sw.now=sw$SweepID, .combine=rbind) %do% {
   print(sw.now)
   wt1 <- wilcox.test(gg.dt[Type=="YRI included" & SweepID==sw.now & !(is.na(pRightTail)),
                            -log10(pRightTail)],
                      gg.dt[Type=="YRI included" & SweepID=="Background"  & !(is.na(pRightTail)),
                            -log10(pRightTail)],
                      paired=F, alternative="greater")
   wt2 <- wilcox.test(gg.dt[Type=="YRI excluded" & SweepID==sw.now & !(is.na(pRightTail)),
                            -log10(pRightTail)],
                      gg.dt[Type=="YRI excluded" & SweepID=="Background"  & !(is.na(pRightTail)),
                            -log10(pRightTail)],
                      paired=F, alternative="greater")
   data.table(SweepID=sw.now, P_YRI=wt1$p.value, P_noYRI=wt2$p.value)

}


##--------------------------------------------------##
# plots
##--------------------------------------------------##

#load("~/Box/Human_adaptation/Datasets/FST.results.OLIVIA.rdata")

#DT <- dt.list$`YRI included`$DT
#ss <- strsplit(DT$SNP, "_")
#DT[, Pos:=as.numeric(sapply(ss, "[", 2))]

# gg.dt <- foreach(k=1:length(dt.list), .combine=rbind) %do% {
#    DT <- data.table(Type=names(dt.list)[k],
#                     dt.list[[k]]$DT[, .(SNP, He, FST,
#                                         OutlierFlag, KEEP,
#                                         pRightTail, qRightTail)])
#
#    ss <- stri_split(DT$SNP, fixed="_")
#    DT[, Chr:=sapply(ss, "[", 1)]
#    DT[, Pos:=as.numeric(sapply(ss, "[", 2))]
#
#    sw.f <- foreach(chr.now=1:22, .combine=rbind) %do% {
#       print(paste("Generating for Chr", chr.now))
#       sw.now <- sw[Chr==chr.now]
#       DT.now <- DT[Chr==chr.now]
#
#       subject <- IRanges(DT.now$Pos, width=rep(1, nrow(DT.now)))
#       query <- IRanges(sw.now$Start, width=sw.now[, End-Start])
#
#       ol <- findOverlaps(query, subject)
#
#       data.table(DT.now[ol@to], sw.now[ol@from, .(SweepID)])
#    }
#    out <- rbind(sw.f,
#                 data.table(DT[!(SNP %in% sw.f$SNP)],
#                            SweepID="Background"))
#
#    out
# }
#
# #SNPs not in sweeps
# mm <- gg.dt[, .(MM=mean(-log10(pRightTail), na.rm=T)),
#             by=c("Type", "SweepID")]
#
# new.ord <- mm[Type=="YRI included"][order(MM), SweepID]
# gg.dt[, SweepID:=factor(SweepID, new.ord)]
# gg.dt[, COL:=factor(ifelse(SweepID=="Background", 1, 0))]
#
# gg.dt[, neglog10p:=-log10(pRightTail)]
# quid <- gg.dt[!is.na(neglog10p), .(Lower=quantile(neglog10p, 0.25),
#                                    Upper=quantile(neglog10p, 0.75),
#                                    Median=median(neglog10p)), by=c("SweepID", "Type")]

# FstDistPlotter(dt.outliers$DF[length(dt.outliers$DF)],
#                DT[KEEP==T, FST],
#                dt.outliers$FSTbar[length(dt.outliers$FSTbar)],
#                binwidth=0.005, titletext=NULL)
#
# DT[KEEP==T, hist(pRightTail, breaks=seq(0, 1, 0.01))]

# xx <- data.table(Population=t(gg.dt[1, 1:13, with=F]),
#                  gg.dt[, .(Mean=sapply(.SD, mean, na.rm=T)),
#                        by="SweepID", .SDcols=grep("^f", colnames(gg.dt))],
#                  gg.dt[, .(serr=sapply(.SD, function(x) sd(x, na.rm=T)/sqrt(.N))),
#                        by="SweepID", .SDcols=grep("^f", colnames(gg.dt))][, .(SE=serr)])
# setnames(xx, c("Population", "SweepID", "Mean", "SE"))
#
# mm <- xx[, .(MM=mean(Mean, na.rm=T)), by="SweepID"]
# new.ord <- mm[order(MM), SweepID]
# xx[, SweepID:=factor(SweepID, new.ord)]
# xx[, COL:=factor(ifelse(SweepID=="Background", 1, 0))]
# xx[, COL2:=factor(ifelse(Population=="YRI", 1, 0))]
#
# ggplot(data=xx[Population!="YRI"],
#        aes(x=SweepID, y=Mean, fill=COL)) +
#    #geom_boxplot() +
#    geom_violin() +
#    geom_point(data=xx[Population=="YRI"], aes(x=SweepID, y=Mean)) +
#    #(yintercept=xx[SweepID=="Background", mean(Mean)], color="red", size=0.25) +
#    #geom_hline(yintercept=xx[SweepID=="Background", median(Mean)], color="gray10", size=0.25) +
#    scale_fill_manual(values=c("lightblue", "white")) +
#    stat_summary(fun.y=median, colour="gray10", geom="point",
#                 shape=18, size=3, show.legend=FALSE) +
#    stat_summary(fun.y=mean, colour="darkred", geom="point",
#                 shape=18, size=3, show.legend=FALSE) +
#    stat_summary(fun.data=mean_se, geom="errorbar", colour="darkred") +
#    theme_bw() +
#    theme(axis.text.x=element_text(angle=45, hjust=1),
#          axis.title.x=element_blank(),
#          legend.position="none")


##--------------------------------------------------##
# Compare corrected and uncorrected FST
##--------------------------------------------------##

# dt.fst[, N.bin:=cut(N1, c(1:20, Inf),
#                     as.character(c(1:19, ">=20")),
#                     right=F, include.lowest=T)]
# dt.fst[, FST.bin:=cut(FST, c(seq(0, 20, 2)/100, 1), include.lowest=T)]
# dt.fst[, He.bin:=cut(He, c(seq(0, 20, 2)/100, 1), include.lowest=T)]
#
# dd <- dt.fst[!(is.na(FST.bin) | is.na(N.bin)), .(N.bin, FST.bin, He.bin,
#                                                  FST, Diff=FST-FST.corr)]
# ggplot(dd, aes(x=FST.bin, y=Diff)) +
#    geom_boxplot() +
#    facet_wrap(~N.bin, ncol=4, drop=TRUE) +
#    #scale_y_log10() +
#    theme_bw() +
#    geom_hline(yintercept=0, color="red", size=0.5) +
#    theme(axis.text.x=element_text(angle=45, hjust=1),
#          legend.position="none")
#
# ggplot(dd, aes(x=N.bin, y=FST)) +
#    geom_boxplot() +
#    #facet_wrap(~N.bin, ncol=4, drop=TRUE) +
#    #scale_y_log10() +
#    theme_bw() +
#    geom_hline(yintercept=0, color="red", size=0.5) +
#    theme(axis.text.x=element_text(angle=45, hjust=1),
#          legend.position="none")
#
# ggplot(dd, aes(x=He.bin, y=Diff)) +
#    geom_boxplot() +
#    facet_wrap(~N.bin, ncol=4, drop=TRUE) +
#    #scale_y_log10() +
#    theme_bw() +
#    geom_hline(yintercept=0, color="red", size=0.5) +
#    theme(axis.text.x=element_text(angle=45, hjust=1),
#          legend.position="none")
#
# #clean up
# rm(dd)
# dt.fst[, N.bin:=NULL]
# dt.fst[, FST.bin:=NULL]
# dt.fst[, He.bin:=NULL]
