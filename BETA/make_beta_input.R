## Cesare de Filippo
## 10-11-2014
## T1 and T2 from DeGiorgio
#modified by Barbara Bitarello on 27.11.2014

library(parallel)
library(SOAR)

source("/home/cesare_filippo/scripts/R_scr/tools/ms_tools.R")
setwd("/mnt/sequencedb/PopGen/barbara/simulations/msms/")



sims.s <- read.ms("Tbs5_f0.5_n100_Sbs0.01.msms.gz",multicore=TRUE, Npop=4,Nchr=c(100,100,100,1)) 
sims.n <- read.ms("neutral_n100.msms.gz",multicore=TRUE, Npop=4,Nchr=c(100,100,100,1)) 
sims.f0.4 <- read.ms("Tbs5_f0.4_n100_Sbs0.01.msms.gz",multicore=TRUE, Npop=4,Nchr=c(100,100,100,1))
sims.f0.3 <- read.ms("Tbs5_f0.3_n100_Sbs0.01.msms.gz",multicore=TRUE, Npop=4,Nchr=c(100,100,100,1))
sims.f0.2 <- read.ms("Tbs5_f0.2_n100_Sbs0.01.msms.gz",multicore=TRUE, Npop=4,Nchr=c(100,100,100,1))
#sims.f0.1 <- read.ms("Tbs5_f0.1_n100_Sbs0.01.msms.gz",multicore=TRUE, Npop=4,Nchr=c(100,100,100,1))

tbs3.sims.s <- read.ms("Tbs3_f0.5_n100_Sbs0.01.msms.gz",multicore=TRUE, Npop=4,Nchr=c(100,100,100,1))
tbs3.sims.n <- read.ms("neutral_n100.msms.gz",multicore=TRUE, Npop=4,Nchr=c(100,100,100,1))
tbs3.sims.f0.4 <- read.ms("Tbs3_f0.4_n100_Sbs0.01.msms.gz",multicore=TRUE, Npop=4,Nchr=c(100,100,100,1))
tbs3.sims.f0.3 <- read.ms("Tbs3_f0.3_n100_Sbs0.01.msms.gz",multicore=TRUE, Npop=4,Nchr=c(100,100,100,1))
tbs3.sims.f0.2 <- read.ms("Tbs3_f0.2_n100_Sbs0.01.msms.gz",multicore=TRUE, Npop=4,Nchr=c(100,100,100,1))


ms2CombinedSNPFile <- function(x, pop="p1",Length=15000) {
  o <- x[grep("p4_1",rownames(x)),]
  ids <- grep(pop, rownames(x))
  bs.pos <- grep("s0.50000",colnames(x))
  a <- apply(x[ids,],2, function(y) sum(y =="0") )
  ## i1 <- which(a == length(ids))
  i1 <- which(a == length(ids) & o == "1")
  a[i1]=0
  ## a <- a[-i1]; o=o[-i1]
  a <- a[a < length(ids)]
  p <- round(as.numeric(sub("s","", names(a)))*Length)
  while(length(unique(p)) != length(p)) {
    i1 <- which(duplicated(p))
    p[i1] <- p[i1-1]+1
  }
#derived
	a<-100-a
	a<-a[a < length(ids)]
  ## ## sample uniformly 
  ## bs.pos <- grep("s0.50000", names(a))
  ## if (bs.pos-1 != length(a) -bs.pos) {
  ##   if(bs.pos-1 > length(a) -bs.pos) {
  ##     i2=(bs.pos-1) - (length(a) -bs.pos)
  ##     i2=1:length(i2);
  ##   } else {

  ##   }
  ##   a <- a[-i2];p=p[-i2]
  ## } 
#  rho <- sapply(2:length(p), function(i) (p[i]-p[i-1])*4*Ne*r)
  out <- list(SNPFile=cbind(position=p,x=a,n=rep(length(ids),length(a))))
} 

## write 10 balancing selection simulations
setwd('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/tmp')


#apply function to all sims
a=mclapply(sims.s, function(x) ms2CombinedSNPFile(x))
a2=mclapply(sims.f0.4, function(x) ms2CombinedSNPFile(x))
a3=mclapply(sims.f0.3, function(x) ms2CombinedSNPFile(x))
a4=mclapply(sims.f0.2, function(x) ms2CombinedSNPFile(x))
#a5=mclapply(sims.f0.1, function(x) ms2CombinedSNPFile(x))

b=n=mclapply(sims.n, function(x) ms2CombinedSNPFile(x))

for( i in 2:2000) {
n[[i]][[1]][,1] = n[[i]][[1]][,1]+1000000*i
} 
#we only need 1000 sim

n=do.call("rbind", lapply(n[1:1000], function(x) x[[1]]))

write.table(n, 'CombinedSNPFile_neu', row.names=F, quote=F, sep="\t")

sapply(1:1000, function(x) sapply(names(b[[x]]), function(y) write.table(b[[x]][[y]], paste0("neu",y,x),row.names=F,col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(a[[x]]), function(y) write.table(a[[x]][[y]], paste0("bs_f0.5_",y,x),row.names=F,col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(a2[[x]]), function(y) write.table(a2[[x]][[y]], paste0("bs_f0.4_",y,x),row.names=F, col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(a3[[x]]), function(y) write.table(a3[[x]][[y]], paste0("bs_f0.3_",y,x),row.names=F, col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(a4[[x]]), function(y) write.table(a4[[x]][[y]], paste0("bs_f0.2_",y,x),row.names=F,col.names=F,quote=F,sep="\t")))
#sapply(1:1000, function(x) sapply(names(a5[[x]]), function(y) write.table(a5[[x]][[y]], paste0("bs_f0.1_",y,x),row.names=F,quote=F,sep="\t")))


#Europe


a_Eu=mclapply(sims.s, function(x) ms2CombinedSNPFile(x, pop='p2'))
a_Eu2=mclapply(sims.f0.4, function(x) ms2CombinedSNPFile(x, pop='p2'))
a_Eu3=mclapply(sims.f0.3, function(x) ms2CombinedSNPFile(x, pop='p2'))
a_Eu4=mclapply(sims.f0.2, function(x) ms2CombinedSNPFile(x, pop='p2'))
#a5=mclapply(sims.f0.1, function(x) ms2CombinedSNPFile(x))

b_Eu=n_Eu=mclapply(sims.n, function(x) ms2CombinedSNPFile(x, pop='p2'))

for( i in 2:2000) {
n_Eu[[i]][[1]][,1] = n_Eu[[i]][[1]][,1]+1000000*i
}
#we only need 1000 sim

n_Eu=do.call("rbind", lapply(n_Eu[1:1000], function(x) x[[1]]))

write.table(n_Eu, 'CombinedSNPFile_neu', row.names=F,col.names=F, quote=F, sep="\t")

sapply(1:1000, function(x) sapply(names(b_Eu[[x]]), function(y) write.table(b_Eu[[x]][[y]], paste0("Eu_neu",y,x),row.names=F,col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(a_Eu[[x]]), function(y) write.table(a_Eu[[x]][[y]], paste0("Eu_bs_f0.5_",y,x),row.names=F,col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(a_Eu2[[x]]), function(y) write.table(a_Eu2[[x]][[y]], paste0("Eu_bs_f0.4_",y,x),row.names=F, col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(a_Eu3[[x]]), function(y) write.table(a_Eu3[[x]][[y]], paste0("Eu_bs_f0.3_",y,x),row.names=F, col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(a_Eu4[[x]]), function(y) write.table(a_Eu4[[x]][[y]], paste0("Eu_bs_f0.2_",y,x),row.names=F,col.names=F,quote=F,sep="\t")))
#sapply(1:1000, function(x) sapply(names(a5[[x]]), function(y) write.table(a5[[x]][[y]], paste0("bs_f0.1_",y,x),row.names=F,quote=F,sep="\t")))

#Asia



a_As=mclapply(sims.s, function(x) ms2CombinedSNPFile(x, pop='p3'))
a_As2=mclapply(sims.f0.4, function(x) ms2CombinedSNPFile(x, pop='p3'))
a_As3=mclapply(sims.f0.3, function(x) ms2CombinedSNPFile(x, pop='p3'))
a_As4=mclapply(sims.f0.2, function(x) ms2CombinedSNPFile(x, pop='p3'))
#a5=mclapply(sims.f0.1, function(x) ms2CombinedSNPFile(x))

b_As=n_As=mclapply(sims.n, function(x) ms2CombinedSNPFile(x, pop='p3'))

for( i in 2:2000) {
n_As[[i]][[1]][,1] = n_As[[i]][[1]][,1]+1000000*i
}
#we only need 1000 sim

n_As=do.call("rbind", lapply(n_As[1:1000], function(x) x[[1]]))

write.table(n_As, 'CombinedSNPFile_neu', row.names=F,col.names=F, quote=F, sep="\t")

sapply(1:1000, function(x) sapply(names(b_As[[x]]), function(y) write.table(b_As[[x]][[y]], paste0("As_neu",y,x),row.names=F,col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(a_As[[x]]), function(y) write.table(a_As[[x]][[y]], paste0("As_bs_f0.5_",y,x),row.names=F,col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(a_As2[[x]]), function(y) write.table(a_As2[[x]][[y]], paste0("As_bs_f0.4_",y,x),row.names=F, col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(a_As3[[x]]), function(y) write.table(a_As3[[x]][[y]], paste0("As_bs_f0.3_",y,x),row.names=F, col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(a_As4[[x]]), function(y) write.table(a_As4[[x]][[y]], paste0("As_bs_f0.2_",y,x),row.names=F,col.names=F,quote=F,sep="\t")))
#sapply(1:1000, function(x) sapply(names(a5[[x]]), function(y) write.table(a5[[x]][[y]], paste0("bs_f0.1_",y,x),row.names=F,quote=F,sep="\t")))



##
setwd('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/BETA/tmp/tbs3/')

#apply function to all sims
tbs3_a=mclapply(tbs3.sims.s, function(x) ms2CombinedSNPFile(x))
tbs3_a2=mclapply(tbs3.sims.f0.4, function(x) ms2CombinedSNPFile(x))
tbs3_a3=mclapply(tbs3.sims.f0.3, function(x) ms2CombinedSNPFile(x))
tbs3_a4=mclapply(tbs3.sims.f0.2, function(x) ms2CombinedSNPFile(x))
#a5=mclapply(sims.f0.1, function(x) ms2CombinedSNPFile(x))

tbs3_b=tbs3_n=mclapply(tbs3.sims.n, function(x) ms2CombinedSNPFile(x))

for( i in 2:2000) {
tbs3_n[[i]][[1]][,1] = tbs3_n[[i]][[1]][,1]+1000000*i
}
#we only need 1000 sim

n=do.call("rbind", lapply(n[1:1000], function(x) x[[1]]))

write.table(n, 'CombinedSNPFile_neu', row.names=F, quote=F, sep="\t")

sapply(1:1000, function(x) sapply(names(tbs3_b[[x]]), function(y) write.table(tbs3_b[[x]][[y]], paste0("neu",y,x),row.names=F,col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(tbs3_a[[x]]), function(y) write.table(tbs3_a[[x]][[y]], paste0("bs_f0.5_",y,x),row.names=F,col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(tbs3_a2[[x]]), function(y) write.table(tbs3_a2[[x]][[y]], paste0("bs_f0.4_",y,x),row.names=F, col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(tbs3_a3[[x]]), function(y) write.table(tbs3_a3[[x]][[y]], paste0("bs_f0.3_",y,x),row.names=F, col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(tbs3_a4[[x]]), function(y) write.table(tbs3_a4[[x]][[y]], paste0("bs_f0.2_",y,x),row.names=F,col.names=F,quote=F,sep="\t")))
#sapply(1:1000, function(x) sapply(names(a5[[x]]), function(y) write.table(a5[[x]][[y]], paste0("bs_f0.1_",y,x),row.names=F,quote=F,sep="\t")))


#Europe
tbs3_a_Eu=mclapply(tbs3.sims.s, function(x) ms2CombinedSNPFile(x, pop='p2'))
tbs3_a2_Eu=mclapply(tbs3.sims.f0.4, function(x) ms2CombinedSNPFile(x, pop='p2'))
tbs3_a3_Eu=mclapply(tbs3.sims.f0.3, function(x) ms2CombinedSNPFile(x, pop='p2'))
tbs3_a4_Eu=mclapply(tbs3.sims.f0.2, function(x) ms2CombinedSNPFile(x, pop='p2'))
#a5=mclapply(sims.f0.1, function(x) ms2CombinedSNPFile(x))

tbs3_b_Eu=tbs3_n_Eu=mclapply(tbs3.sims.n, function(x) ms2CombinedSNPFile(x, pop='p2'))

for( i in 2:2000) {
tbs3_n_Eu[[i]][[1]][,1] = tbs3_n_Eu[[i]][[1]][,1]+1000000*i
}
#we only need 1000 sim


sapply(1:1000, function(x) sapply(names(tbs3_b_Eu[[x]]), function(y) write.table(tbs3_b_Eu[[x]][[y]], paste0("Eu_neu",y,x),row.names=F,col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(tbs3_a_Eu[[x]]), function(y) write.table(tbs3_a_Eu[[x]][[y]], paste0("Eu_bs_f0.5_",y,x),row.names=F,col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(tbs3_a2_Eu[[x]]), function(y) write.table(tbs3_a2_Eu[[x]][[y]], paste0("Eu_bs_f0.4_",y,x),row.names=F, col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(tbs3_a3_Eu[[x]]), function(y) write.table(tbs3_a3_Eu[[x]][[y]], paste0("Eu_bs_f0.3_",y,x),row.names=F, col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(tbs3_a4_Eu[[x]]), function(y) write.table(tbs3_a4_Eu[[x]][[y]], paste0("Eu_bs_f0.2_",y,x),row.names=F,col.names=F,quote=F,sep="\t")))
#sapply(1:1000, function(x) sapply(names(a5[[x]]), function(y) write.table(a5[[x]][[y]], paste0("bs_f0.1_",y,x),row.names=F,quote=F,sep="\t")))

#ASIA


tbs3_a_As=mclapply(tbs3.sims.s, function(x) ms2CombinedSNPFile(x, pop='p3'))
tbs3_a2_As=mclapply(tbs3.sims.f0.4, function(x) ms2CombinedSNPFile(x, pop='p3'))
tbs3_a3_As=mclapply(tbs3.sims.f0.3, function(x) ms2CombinedSNPFile(x, pop='p3'))
tbs3_a4_As=mclapply(tbs3.sims.f0.2, function(x) ms2CombinedSNPFile(x, pop='p3'))
#a5=mclapply(sims.f0.1, function(x) ms2CombinedSNPFile(x))

tbs3_b_As=tbs3_n_As=mclapply(tbs3.sims.n, function(x) ms2CombinedSNPFile(x, pop='p3'))

for( i in 2:2000) {
tbs3_n_As[[i]][[1]][,1] = tbs3_n_As[[i]][[1]][,1]+1000000*i
}
#we only need 1000 sim


sapply(1:1000, function(x) sapply(names(tbs3_b_As[[x]]), function(y) write.table(tbs3_b_As[[x]][[y]], paste0("As_neu",y,x),row.names=F,col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(tbs3_a_As[[x]]), function(y) write.table(tbs3_a_As[[x]][[y]], paste0("As_bs_f0.5_",y,x),row.names=F,col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(tbs3_a2_As[[x]]), function(y) write.table(tbs3_a2_As[[x]][[y]], paste0("As_bs_f0.4_",y,x),row.names=F, col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(tbs3_a3_As[[x]]), function(y) write.table(tbs3_a3_As[[x]][[y]], paste0("As_bs_f0.3_",y,x),row.names=F, col.names=F,quote=F,sep="\t")))
sapply(1:1000, function(x) sapply(names(tbs3_a4_As[[x]]), function(y) write.table(tbs3_a4_As[[x]][[y]], paste0("As_bs_f0.2_",y,x),row.names=F,col.names=F,quote=F,sep="\t")))
#sapply(1:1000, function(x) sapply(names(a5[[x]]), function(y) write.table(a5[[x]][[y]], paste0("bs_f0.1_",y,x),row.names=F,quote=F,sep="\t")))



