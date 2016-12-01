###################################
#
#	Barbara Bitarello
#
#	Last modified: 01/12/2016
#
##################################


library(parallel)
library(SOAR)
library(ggplot2)
library(vioplot)
library(dplyr)
library(reshape)
library(reshape2)
Sys.setenv(R_LOCAL_CACHE="estsession")
library(data.table)
library(qqman)

pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")


library(grDevices)

a<-c('[0-1]', ']1-2]', ']2-3]', ']3-4]', ']4-5]', '>5')

for.plot<-c(list.SCAN[[2]]$PtoD[which(list.SCAN[[2]]$PtoD<=5)], rep(6, sum(list.SCAN[[2]]$PtoD>5)))

for.plot2<-c(Union.CANDf0.5_0.4_0.3[[2]]$PtoD[which(Union.CANDf0.5_0.4_0.3[[2]]$PtoD<=5)], rep(6, sum(Union.CANDf0.5_0.4_0.3[[2]]$PtoD>5)))
h<-hist(for.plot, nclass=6,plot=F)
hist(for.plot2,breaks=seq(0,6,1), plot=F)-> h1

pdf('figures/PtoD.for.paper.LWK.pdf')
adjustcolor('darkgray', alpha.f=0.50)-> coor_transparent
adjustcolor('darkorchid3', alpha.f=0.5)-> coor_transparent2
barplot((h$counts/length(for.plot)),col=coor_transparent,space=1, border=F, ylim=c(0, 0.9), axes=F)->bp
barplot((h1$counts/length(for.plot2)),col=coor_transparent2,border=F,space=1,  add=T, axes=F)->bp2
axis(1,at=c(bp),labels=a, cex.axis=1.2, las=2)
axis(2, cex.axis=1.2)
title(ylab="Relative Frequency",xlab="", main="SNP/(FD+1)", cex.lab=1.2)
legend('topright', c("Background", "Significant Windows"), col=c(coor_transparent, coor_transparent2), pch=c(22,22), fill=c(coor_transparent, coor_transparent2),bty="n")
dev.off()


for.plot.GBR<-c(list.SCAN[[6]]$PtoD[which(list.SCAN[[6]]$PtoD<=5)], rep(5, sum(list.SCAN[[6]]$PtoD>5)))

for.plot2.GBR<-c(Union.CANDf0.5_0.4_0.3[[6]]$PtoD[which(Union.CANDf0.5_0.4_0.3[[6]]$PtoD<=5)], rep(6, sum(Union.CANDf0.5_0.4_0.3[[6]]$PtoD>5)))
h1GBR<-hist(for.plot2.GBR, plot=F,breaks=seq(0,6, 1))

hGBR<-hist(for.plot.GBR, plot=F, breaks=seq(0,6,1))

pdf('figures/PtoD.for.paper.GBR.pdf')
adjustcolor('darkgray', alpha.f=0.50)-> coL1
adjustcolor('darkolivegreen', alpha.f=0.5)-> coL2
barplot((hGBR$counts/length(for.plot)),col=coL1,space=1, border=F, ylim=c(0, 0.9), axes=F)->bp

barplot((h1GBR$counts/length(for.plot2)),col=coL2,border=F,space=1,  add=T, axes=F)->bp2

axis(1,at=c(bp),labels=a, cex.axis=1.2, las=2)
axis(2, cex.axis=1.2)
title(ylab="Relative Frequency",xlab="", main="SNP/(FD+1)", cex.lab=1.2)

legend('topright', c("Background", "Significant Windows"), col=c(coL1, coL2), pch=c(22,22), fill=c(coL1, coL2),bty="n")

dev.off()

