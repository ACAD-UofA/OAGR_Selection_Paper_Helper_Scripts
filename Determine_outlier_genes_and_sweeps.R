library(gdata)
library(data.table)
library(matrixStats)
library(foreach)
library(dplyr)
library(stringi)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(IRanges)

#set working directory
work.dir <- "~/Dropbox/Selection_test/"
setwd(work.dir)

#path to raw files (genes + pathway info)
raw.path <- file.path(work.dir, "Raw_files")

#path to sweepfinder files
sf.path <- file.path(work.dir, "SweepFinder")

# load PolyLink libraries
scripts.path <- file.path(work.dir, "Scripts")
source(file.path(scripts.path, "PolyLink.R")) 

#make output directory
out.path <- file.path(work.dir, "GeneScores")
dir.create(out.path)

##----------------------------------------
#read in gene data
##----------------------------------------

#read in list of genes
prots <- fread(file.path(raw.path, "Entrez_Gene_IDs_positions_HG19.txt"))

#add upper and lower flanking region boundaries (set to 0 if exact gene boundaries desired)
flank <- 50000
prots[, Lower:=ifelse((Start-flank)<0, 0, Start-flank)]
prots[, Upper:=End+flank]
n.prots <- nrow(prots)

##----------------------------------------
#collate SweepFinder results by gene
##----------------------------------------

sf.f <- dir(sf.path, no.. = T)

#read in one chromosome at a time for each population
for(chr in 1:22){ # loop over each chromosome
   print(paste("reading", d, "data, chr", chr))
   sf <- foreach(f=sf.f, .combine=rbind) %do% { # loop over all pops in given classification
      print(paste(chr, f))
      ll <- list.files(file.path(sf.path, f, "out"))
      f.now <- ll[grep(paste0("_", chr, "\\."), ll)]
      if(length(f.now)==0){
         print(paste("Missing", f, "chr", chr))
         dd <- NULL
      }else{
         ff <- fread(file.path(sf.path, f, "out", f.now))
         dd <- data.table(Population=f, ff[, 1:2, with=F])
      }
      dd
   }

   #get chromosome specific genes
   prot.now <- prots[Chrom==chr]

   #find intersecting scores for each gene (using IRanges)
   print("Calculating optimal CLR values")
   subject <- IRanges(sf$location, width=rep(1, nrow(sf)))
   query <- IRanges(prot.now$Lower, width=prot.now[, Upper-Lower])

   ol <- findOverlaps(query, subject)

   score.dt <- data.table(sf[ol@to], 
                          prot.now[ol@from, .(Chr=Chrom, GeneID, 
                                              Start, End, Biotype)])

   out <- score.dt[, .(LR.gene.max=max(LR, na.rm=T)),
                    by=c("Population", "GeneID", "Chr",
                         "Start", "End", "Biotype")]

   #include missing genes (if any)
   unq.genes <- prot.now[, .(GeneID, Chr=Chrom, Start, End, Biotype)]
   missing.genes <- foreach(pop=out[, unique(Population)], .combine=rbind) %do% {
      missing <- out[Population==pop, setdiff(unq.genes$GeneID, GeneID)]
      if(length(missing>0)) data.table(Population=pop, 
                                       unq.genes[GeneID %in% missing], 
                                       LR.gene.max=NA)
   }
   out <- rbind(out, missing.genes)[order(Population, GeneID)]

   print(paste("writing", d, "chr", chr))
   fwrite(out, file.path(out.path, "GeneScores_max_LR.txt"),
          na="NA", sep="\t", quote=F, append=T)
}


##----------------------------------------
#add log transformed scores
##----------------------------------------

print("Including log transformed, gene length corrected scores")
all.sf <- fread(file.path(out.path, "GeneScores_max_LR.txt"))

#Log transform scores
all.sf[, log10.LR.gene.max:=log10(LR.gene.max+0.01)]

#adjust score for gene length (uses function from PolyLink -- SLOW)
minbs <- 2000; maxbs <- 2500 #defaults
all.sf[, GeneLength:=End-Start]

unq.genes <- all.sf[!duplicated(GeneID), .(GeneID, GeneLength)]

bins <- AssignBins(unq.genes, fld="GeneLength", 
                   min.bin.size=minbs, max.bin.size=maxbs)

out <- merge(all.sf, bins[, .(GeneID, objBin)], by="GeneID")


#modified function to allow for NAs (adapted from PolyLink function, allowing for missing data)
modz<- function(x){
   y<-x-median(x, na.rm=T)
   return(0.6745*y/median(abs(y), na.rm=T))
}

out[, objScore:=modz(log10.LR.gene.max), by=c("Population", "objBin")]

print(paste("writing", d)) 
fwrite(out, file.path(out.path, "GeneScores_max_LR_GeneLengthCorrected.txt"),
       na="NA", sep="\t", quote=F)

##-----------------------------------------------##
# calculate p and q values
##-----------------------------------------------##

dt <- fread(file.path(out.path, "GeneScores_max_LR_GeneLengthCorrected.txt"))

dt[, objScore.2tail.p:=2*(1-pnorm(abs(objScore), mean=0, sd=1))]
dt[objScore>0, objScore.gt0.fdr.2tail.q:=qvalue(objScore.2tail.p)$qvalues, 
   by="Population"]


##-----------------------------------------------##
# add effective population sizes
##-----------------------------------------------##

#read in one chromosome at a time for each population
sfs <- foreach(chr=1:22, .combine=rbind) %do% { # loop over each chromosome
   print(paste("reading", d, "data, chr", chr))
   sf <- foreach(f=sf.f, .combine=rbind) %do% { # loop over all pops in given classification
      print(paste(chr, f))
      ll <- list.files(file.path(sf.path, f, "sfs"))
      f.now <- ll[grep(paste0("_", chr, "\\."), ll)]
      if(length(f.now)==0){
         print(paste("Missing", f, "chr", chr))
         dd <- NULL
      }else{
         ff <- fread(file.path(sf.path, f, "sfs", f.now))
         dd <- data.table(Population=f, ff)
      }
      dd
   }
}

samp.info <- sfs[, .(samp.N=max(n), eff.N=mean(n), 
                     var.eff.N=var(n)), by="Population"]

fwrite(samp.info, col.names=T, row.names=F, quote=F, sep='\t',
       file=file.path(out.path, "effective.N.info.txt"))


##-----------------------------------------------##
# diagnostic plots: Z and p value distributions
##-----------------------------------------------##

#add sample size information
dt.m <- merge(dt, samp.info, by="Population")

#add new label
dt.m[, Pop.name:=paste0(Population, " [", round(eff.N, 2), "]")]
ord.eff.n <- dt.m[!duplicated(Population)][order(eff.N), 
                                           paste0(Population, " [", round(eff.N, 2), "]")]
dt.m[, Pop.name:=factor(Pop.name, levels=ord.eff.n)]

plot.path <- file.path(work.dir, "Plots")
dir.create(plot.path)

#number of columns in plots (alter to preferences)
n <- dt.m[, length(unique(Population))]

#Z score precorrection
ggplot(dt.m, aes(GeneLength, log10.LR.gene.max)) +
   geom_hex() +
   geom_smooth(method=lm) + scale_x_log10() + theme_bw() +
   scale_fill_gradient(name="count", trans="log10") +
   ylab("Uncorrected Scores") +
   facet_wrap(~Pop.name, ncol=n) +
   xlab("Population") +
   theme_bw() +
   theme(axis.text.x=element_text(size=6, angle=45, hjust=1),
         axis.text.y=element_text(size=6),
         strip.text=element_text(size=6),
         panel.border=element_blank(),
         panel.spacing.x=unit(0.025, "lines"),
         panel.spacing.y=unit(0.01, "lines"),
         legend.position="bottom",
         legend.direction="horizontal",
         legend.title=element_text(size=6),
         legend.text=element_text(size=5.5))

ggsave(filename=file.path(plot.path, "ObjStat_vs_GeneLength_precorrection.pdf"),
       width=8, height=6, dpi=300)

#Z score postcorrection
ggplot(dt.m, aes(GeneLength, objScore)) +
   geom_hex() +
   geom_smooth(method=lm) + scale_x_log10() + theme_bw() +
   scale_fill_gradient(name="count", trans="log10") +
   facet_wrap(~Pop.name, ncol=n) +
   xlab("Population") + ylab("Corrected Scores") +
   theme_bw() +
      theme(axis.text.x=element_text(size=6, angle=45, hjust=1),
         axis.text.y=element_text(size=6),
         strip.text=element_text(size=6),
         panel.border=element_blank(),
         panel.spacing.x=unit(0.025, "lines"),
         panel.spacing.y=unit(0.01, "lines"),
         legend.position="botton",
         legend.direction="horizontal",
         legend.title=element_text(size=6),
         legend.text=element_text(size=5.5))

ggsave(filename=file.path(plot.path, "ObjStat_vs_GeneLength_postcorrection.pdf"),
       width=8, height=6, dpi=300)

#Z scores
ggplot(dt.m, aes(x=objScore)) +
   geom_histogram(aes(y =..density..), binwidth=0.5,
                  colour="black", fill="white", size=0.25) +
   facet_wrap(~Pop.name, ncol=n) +
   stat_function(fun=dnorm, args=list(mean=0, sd=1), col="red") +
   xlab("Gene Z score") + ylab("Density") +
   theme_bw() +
   theme(axis.text.x=element_text(size=6, hjust=1),
         axis.text.y=element_text(size=6),
         strip.text=element_text(size=8),
         panel.border=element_blank(),
         panel.spacing.x=unit(0.025, "lines"),
         panel.spacing.y=unit(0.01, "lines"))

ggsave(filename=file.path(plot.path, "Obj.scores.distribution.pdf"),
       width=8, height=6, dpi=300)

#p-values
ggplot(dt.m[objScore>0], aes(x=objScore.2tail.p)) +
   geom_histogram(aes(y =..density..),
                  breaks=seq(0,1,0.01)) +
   facet_wrap(~Pop.name, ncol=n) +
   xlab("p value") + ylab("Density") +
   theme_bw() +
   theme(axis.text.x=element_text(size=6, hjust=1),
         axis.text.y=element_text(size=6),
         strip.text=element_text(size=8),
         panel.border=element_blank(),
         panel.spacing.x=unit(0.025, "lines"),
         panel.spacing.y=unit(0.01, "lines"))

ggsave(filename=file.path(plot.path, "p.values.objScoreGT0.2tailed.distribution.pdf"),
       width=8, height=6, dpi=300)



##-----------------------------------------------##
# determine outliers and bin into sweeps
##-----------------------------------------------##

#do for all populations and excluding small populations
small.pops <- dt.m[floor(eff.N) <= 10, unique(Population)]
dt.m[, Class:="Ancient"]
dt.m[Population %in% small.pops, Class:="Small"]


dt.sweep <- foreach(pops=c("AllPops", "NoSmall")) %do% {
   if(pops=="AllPops") pops.now <- dt.m[, unique(Population)]
   if(pops=="NoSmall") pops.now <- dt.m[!(Population %in% small.pops),
                                              unique(Population)]

   foreach(q.cutoff=c(0.1, 0.05, 0.01), .combine=rbind) %do% {
      #determine outlier genes
      dt.outliers <- dt.m[Population %in% pops.now &
                          objScore.gt0.fdr.2tail.q<q.cutoff][order(Population, Chr, Start)]

      #combine genes less than xx bp distant (check sensitivity of this variable)
      foreach(sweep.join=c(250000, 500000, 1000000), .combine=rbind) %do% {
         print(paste0("Running for ", pops, " q<=", q.cutoff,
                      " & join genes <", sweep.join, " apart"))
         print("Combining sweeps step 1...")
         sweep.comb <- foreach(pop=dt.m[, unique(Population)], .combine=rbind) %do% {
            cc=0
            foreach(chr=dt.m[Population==pop, unique(Chr)], .combine=rbind) %do% {
               dt.out.now <- dt.outliers[Population==pop & Chr==chr]
               if(nrow(dt.out.now)>0){
                  out <- dt.out.now[, .(Population, Chr, Start, End,
                                        SweepBin=cc+cumsum(c(TRUE, !diff((Start+End)/2)<=sweep.join)))]
                  cc=max(out$SweepBin, na.rm=T)
               }else{
                  out <- NULL
               }
               out
            }
         }

         dt.outliers[, Sweep.1M.bin:=sweep.comb[order(Population, Chr, Start), SweepBin]]

         #determine overlapping sweep windows across different populations
         print("Combining sweeps step 2...")

         sweep.width <- dt.outliers[, .(Chr=unique(Chr), Start=min(Start), End=max(End)),
                                    by=c("Population", "Sweep.1M.bin")]

         unq.chr <- sweep.width[, sort(unique(Chr))]
         cc=1 # chromosome number
         for(chr.now in unq.chr){
            sw.n <- sweep.width[Chr==chr.now]
            subject <- IRanges(sw.n$Start, width=sw.n[, End-Start])
            query <- IRanges(sw.n$Start, width=sw.n[, End-Start])

            ol <- findOverlaps(query, subject)

            unq.from <- unique(ol@from)
            while(length(unq.from)>0){
               i=unq.from[1]
               from.to <- ol[ol@from==i]@to
               l=1
               while(l>0){
                  from.to.2 <- unique(ol[ol@from %in% from.to]@to)
                  l <- length(setdiff(from.to.2, from.to))
                  from.to <- from.to.2
               }
               o2 <- sw.n[from.to]
               for(j in 1:nrow(o2)){
                  dt.outliers[Population==o2[j, Population] &
                              Sweep.1M.bin==o2[j, Sweep.1M.bin], Sweep.1M.bin.comb:=cc]
               }
               cc=cc+1
               unq.from <- setdiff(unq.from, from.to)
            }
         }
         #add sweep label
         dt.outliers[, SweepID:=paste0(Chr, ":", round(min(Start)/1000000, 1),
                                       "-", max(round(max(End)/1000000, 1),
                                                round(min(Start)/1000000, 1)+0.1)), 
                     by="Sweep.1M.bin.comb"]
         
         data.table(Join.width=sweep.join, Cutoff=q.cutoff, dt.outliers)
      }
   }
}
names(dt.sweep) <- c("AllPops", "NoSmall")
