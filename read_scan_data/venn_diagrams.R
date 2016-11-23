##########################################
#	Bárbara D Bitarello
#
#	Last modified: 21.11.2-16
#########################################


library(parallel)  #parallelize functions
library(SOAR)  #store objects from workspace
library(ggplot2)  #pretty plots
library(plyr)  #big data frames
library(dplyr)
library(qqman)  #manhattan plot
library(VennDiagram)  
require(mgcv)
library(data.table)


Sys.setenv(R_LOCAL_CACHE="estsession")  #this is for 'SOAR'
pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")



#windows

venn.diagram(list(LWK=Union.CANDf0.5_0.4_0.3[[2]]$Win.ID, YRI=Union.CANDf0.5_0.4_0.3[[3]]$Win.ID, GBR=Union.CANDf0.5_0.4_0.3[[6]]$Win.ID, TSI=Union.CANDf0.5_0.4_0.3[[7]]$Win.ID), fill=c("slateblue4", "slateblue", "darkolivegreen4", "darkolivegreen3"), alpha = c(0.5,0.5,0.5,0.5), cex=2, lty=1, fontfamily=3, filename="figures/VENN.UNION.CAND.tiff",imagetype="png")


venn.diagram(list(LWK=Union.top0.5_0.4_0.3[[2]]$Win.ID, YRI=Union.top0.5_0.4_0.3[[3]]$Win.ID, GBR=Union.top0.5_0.4_0.3[[6]]$Win.ID, TSI=Union.top0.5_0.4_0.3[[7]]$Win.ID), fill=c("slateblue4", "slateblue", "darkolivegreen4", "darkolivegreen3"), alpha = c(0.5,0.5,0.5,0.5), cex=2, lty=1, fontfamily=3, filename="figures/VENN.UNION.top829.tiff",imagetype="png")




#genes

read.table('bedfiles/prot.cod.Union.CAND.0.5_0.4_0.3_LWK.bed')-> LWK.cand.genes
read.table('bedfiles/prot.cod.Union.CAND.0.5_0.4_0.3_YRI.bed')-> YRI.cand.genes
read.table('bedfiles/prot.cod.Union.CAND.0.5_0.4_0.3_GBR.bed')-> GBR.cand.genes
read.table('bedfiles/prot.cod.Union.CAND.0.5_0.4_0.3_TSI.bed')-> TSI.cand.genes



venn.diagram(list(LWK=LWK.cand.genes$V1, YRI=YRI.cand.genes$V1, GBR=GBR.cand.genes$V1, TSI=TSI.cand.genes$V1), fill=c("slateblue4", "slateblue", "darkolivegreen4", "darkolivegreen3"), alpha = c(0.5,0.5,0.5,0.5), cex=2, lty=1, fontfamily=3, filename="figures/VENN.UNION.CAND.genes.tiff",imagetype="png")

read.table('bedfiles/prot.cod.Union.top829.0.5_0.4_0.3_LWK.bed')-> LWK.top.genes
read.table('bedfiles/prot.cod.Union.top829.0.5_0.4_0.3_YRI.bed')-> YRI.top.genes
read.table('bedfiles/prot.cod.Union.top829.0.5_0.4_0.3_GBR.bed')-> GBR.top.genes
read.table('bedfiles/prot.cod.Union.top829.0.5_0.4_0.3_TSI.bed')-> TSI.top.genes



venn.diagram(list(LWK=LWK.top.genes$V1, YRI=YRI.top.genes$V1, GBR=GBR.top.genes$V1, TSI=TSI.top.genes$V1), fill=c("slateblue4", "slateblue", "darkolivegreen4", "darkolivegreen3"), alpha = c(0.5,0.5,0.5,0.5), cex=2, lty=1, fontfamily=3, filename="figures/VENN.UNION.top829.genes.tiff",imagetype="png")

#now each tf separately

venn.diagram(list(LWK=CANDf0.5[[2]]$Win.ID, YRI=CANDf0.5[[3]]$Win.ID, GBR=CANDf0.5[[6]]$Win.ID, TSI=CANDf0.5[[7]]$Win.ID), fill=c("slateblue4", "slateblue", "darkolivegreen4", "darkolivegreen3"), alpha = c(0.5,0.5,0.5,0.5), cex=2, lty=1, fontfamily=3, filename="figures/VENN.CANDf0.5.windows.tiff",imagetype="png")


venn.diagram(list(LWK=CANDf0.4[[2]]$Win.ID, YRI=CANDf0.4[[3]]$Win.ID, GBR=CANDf0.4[[6]]$Win.ID, TSI=CANDf0.4[[7]]$Win.ID), fill=c("slateblue4", "slateblue", "darkolivegreen4", "darkolivegreen3"), alpha = c(0.5,0.5,0.5,0.5), cex=2, lty=1, fontfamily=3, filename="figures/VENN.CANDf0.4.windows.tiff",imagetype="png")


venn.diagram(list(LWK=CANDf0.3[[2]]$Win.ID, YRI=CANDf0.3[[3]]$Win.ID, GBR=CANDf0.3[[6]]$Win.ID, TSI=CANDf0.3[[7]]$Win.ID), fill=c("slateblue4", "slateblue", "darkolivegreen4", "darkolivegreen3"), alpha = c(0.5,0.5,0.5,0.5), cex=2, lty=1, fontfamily=3, filename="figures/VENN.CANDf0.3.windows.tiff",imagetype="png")



venn.diagram(list(LWK=top829f0.5[[2]]$Win.ID, YRI=top829f0.5[[3]]$Win.ID, GBR=top829f0.5[[6]]$Win.ID, TSI=top829f0.5[[7]]$Win.ID), fill=c("slateblue4", "slateblue", "darkolivegreen4", "darkolivegreen3"), alpha = c(0.5,0.5,0.5,0.5), cex=2, lty=1, fontfamily=3, filename="figures/VENN.top829f0.5.windows.tiff",imagetype="png")


venn.diagram(list(LWK=top829f0.4[[2]]$Win.ID, YRI=top829f0.4[[3]]$Win.ID, GBR=top829f0.4[[6]]$Win.ID, TSI=top829f0.4[[7]]$Win.ID), fill=c("slateblue4", "slateblue", "darkolivegreen4", "darkolivegreen3"), alpha = c(0.5,0.5,0.5,0.5), cex=2, lty=1, fontfamily=3, filename="figures/VENN.top829f0.4.windows.tiff",imagetype="png")


venn.diagram(list(LWK=top829f0.3[[2]]$Win.ID, YRI=top829f0.3[[3]]$Win.ID, GBR=top829f0.3[[6]]$Win.ID, TSI=top829f0.3[[7]]$Win.ID), fill=c("slateblue4", "slateblue", "darkolivegreen4", "darkolivegreen3"), alpha = c(0.5,0.5,0.5,0.5), cex=2, lty=1, fontfamily=3, filename="figures/VENN.top829f0.3.windows.tiff",imagetype="png")



#genes

read.table('bedfiles/prot.cod.CANDf0.5.LWK.bed')-> LWK.f0.5.genes
read.table('bedfiles/prot.cod.CANDf0.5.YRI.bed')-> YRI.f0.5.genes
read.table('bedfiles/prot.cod.CANDf0.5.GBR.bed')-> GBR.f0.5.genes
read.table('bedfiles/prot.cod.CANDf0.5.TSI.bed')-> TSI.f0.5.genes

read.table('bedfiles/prot.cod.CANDf0.4.LWK.bed')-> LWK.f0.4.genes
read.table('bedfiles/prot.cod.CANDf0.4.YRI.bed')-> YRI.f0.4.genes
read.table('bedfiles/prot.cod.CANDf0.4.GBR.bed')-> GBR.f0.4.genes
read.table('bedfiles/prot.cod.CANDf0.4.TSI.bed')-> TSI.f0.4.genes

read.table('bedfiles/prot.cod.CANDf0.3.LWK.bed')-> LWK.f0.3.genes
read.table('bedfiles/prot.cod.CANDf0.3.YRI.bed')-> YRI.f0.3.genes
read.table('bedfiles/prot.cod.CANDf0.3.GBR.bed')-> GBR.f0.3.genes
read.table('bedfiles/prot.cod.CANDf0.3.TSI.bed')-> TSI.f0.3.genes

read.table('bedfiles/prot.cod.top829f0.5.LWK.bed')-> LWK.f0.5.top829.genes
read.table('bedfiles/prot.cod.top829f0.5.YRI.bed')-> YRI.f0.5.top829.genes
read.table('bedfiles/prot.cod.top829f0.5.GBR.bed')-> GBR.f0.5.top829.genes
read.table('bedfiles/prot.cod.top829f0.5.TSI.bed')-> TSI.f0.5.top829.genes

read.table('bedfiles/prot.cod.top829f0.4.LWK.bed')-> LWK.f0.4.top829.genes
read.table('bedfiles/prot.cod.top829f0.4.YRI.bed')-> YRI.f0.4.top829.genes
read.table('bedfiles/prot.cod.top829f0.4.GBR.bed')-> GBR.f0.4.top829.genes
read.table('bedfiles/prot.cod.top829f0.4.TSI.bed')-> TSI.f0.4.top829.genes

read.table('bedfiles/prot.cod.top829f0.3.LWK.bed')-> LWK.f0.3.top829.genes
read.table('bedfiles/prot.cod.top829f0.3.YRI.bed')-> YRI.f0.3.top829.genes
read.table('bedfiles/prot.cod.top829f0.3.GBR.bed')-> GBR.f0.3.top829.genes
read.table('bedfiles/prot.cod.top829f0.3.TSI.bed')-> TSI.f0.3.top829.genes



venn.diagram(list(LWK=LWK.f0.5.genes$V1, YRI=YRI.f0.5.genes$V1, GBR=GBR.f0.5.genes$V1, TSI=TSI.f0.5.genes$V1), fill=c("slateblue4", "slateblue", "darkolivegreen4", "darkolivegreen3"), alpha = c(0.5,0.5,0.5,0.5), cex=2, lty=1, fontfamily=3, filename="figures/VENN.CAND.f0.5.genes.tiff",imagetype="png")


venn.diagram(list(LWK=LWK.f0.4.genes$V1, YRI=YRI.f0.4.genes$V1, GBR=GBR.f0.4.genes$V1, TSI=TSI.f0.4.genes$V1), fill=c("slateblue4", "slateblue", "darkolivegreen4", "darkolivegreen3"), alpha = c(0.5,0.5,0.5,0.5), cex=2, lty=1, fontfamily=3, filename="figures/VENN.CAND.f0.4.genes.tiff",imagetype="png")


venn.diagram(list(LWK=LWK.f0.3.genes$V1, YRI=YRI.f0.3.genes$V1, GBR=GBR.f0.3.genes$V1, TSI=TSI.f0.3.genes$V1), fill=c("slateblue4", "slateblue", "darkolivegreen4", "darkolivegreen3"), alpha = c(0.5,0.5,0.5,0.5), cex=2, lty=1, fontfamily=3, filename="figures/VENN.CAND.f0.3.genes.tiff",imagetype="png")

venn.diagram(list(LWK=LWK.f0.5.top829.genes$V1, YRI=YRI.f0.5.top829.genes$V1, GBR=GBR.f0.5.top829.genes$V1, TSI=TSI.f0.5.top829.genes$V1), fill=c("slateblue4", "slateblue", "darkolivegreen4", "darkolivegreen3"), alpha = c(0.5,0.5,0.5,0.5), cex=2, lty=1, fontfamily=3, filename="figures/VENN.top829.f0.5.genes.tiff",imagetype="png")

venn.diagram(list(LWK=LWK.f0.4.top829.genes$V1, YRI=YRI.f0.4.top829.genes$V1, GBR=GBR.f0.4.top829.genes$V1, TSI=TSI.f0.4.top829.genes$V1), fill=c("slateblue4", "slateblue", "darkolivegreen4", "darkolivegreen3"), alpha = c(0.5,0.5,0.5,0.5), cex=2, lty=1, fontfamily=3, filename="figures/VENN.top829.f0.4.genes.tiff",imagetype="png")


venn.diagram(list(LWK=LWK.f0.3.top829.genes$V1, YRI=YRI.f0.3.top829.genes$V1, GBR=GBR.f0.3.top829.genes$V1, TSI=TSI.f0.3.top829.genes$V1), fill=c("slateblue4", "slateblue", "darkolivegreen4", "darkolivegreen3"), alpha = c(0.5,0.5,0.5,0.5), cex=2, lty=1, fontfamily=3, filename="figures/VENN.top829.f0.3.genes.tiff",imagetype="png")



