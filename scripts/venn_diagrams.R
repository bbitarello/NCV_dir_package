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


#calculating overlaps without plotting venn diagrmas

unique(sort(c(intersect(Union.CANDf0.5_0.4_0.3[[6]]$Win.ID,Union.CANDf0.5_0.4_0.3[[2]]$Win.ID),  intersect(Union.CANDf0.5_0.4_0.3[[6]]$Win.ID,Union.CANDf0.5_0.4_0.3[[3]]$Win.ID), intersect(Union.CANDf0.5_0.4_0.3[[6]]$Win.ID,Union.CANDf0.5_0.4_0.3[[7]]$Win.ID)))) #7637/9521=80.21%


unique(sort(c(intersect(Union.CANDf0.5_0.4_0.3[[7]]$Win.ID,Union.CANDf0.5_0.4_0.3[[2]]$Win.ID),  intersect(Union.CANDf0.5_0.4_0.3[[7]]$Win.ID,Union.CANDf0.5_0.4_0.3[[3]]$Win.ID), intersect(Union.CANDf0.5_0.4_0.3[[7]]$Win.ID,Union.CANDf0.5_0.4_0.3[[6]]$Win.ID)))) #7612/9282=82.01%


unique(sort(c(intersect(Union.CANDf0.5_0.4_0.3[[3]]$Win.ID,Union.CANDf0.5_0.4_0.3[[2]]$Win.ID),  intersect(Union.CANDf0.5_0.4_0.3[[3]]$Win.ID,Union.CANDf0.5_0.4_0.3[[6]]$Win.ID), intersect(Union.CANDf0.5_0.4_0.3[[3]]$Win.ID,Union.CANDf0.5_0.4_0.3[[7]]$Win.ID)))) #8302/10997=75.49%

unique(sort(c(intersect(Union.CANDf0.5_0.4_0.3[[2]]$Win.ID,Union.CANDf0.5_0.4_0.3[[3]]$Win.ID),  intersect(Union.CANDf0.5_0.4_0.3[[2]]$Win.ID,Union.CANDf0.5_0.4_0.3[[6]]$Win.ID), intersect(Union.CANDf0.5_0.4_0.3[[2]]$Win.ID,Union.CANDf0.5_0.4_0.3[[7]]$Win.ID)))) #7984/10072=79.27%

unique(sort(c(intersect(Union.CANDf0.5_0.4_0.3[[6]]$Win.ID,Union.CANDf0.5_0.4_0.3[[2]]$Win.ID),  intersect(Union.CANDf0.5_0.4_0.3[[6]]$Win.ID,Union.CANDf0.5_0.4_0.3[[3]]$Win.ID), intersect(Union.CANDf0.5_0.4_0.3[[6]]$Win.ID,Union.CANDf0.5_0.4_0.3[[7]]$Win.ID)))) #7637/9521=80.21%


#within continents

unique(sort(c(intersect(Union.CANDf0.5_0.4_0.3[[3]]$Win.ID,Union.CANDf0.5_0.4_0.3[[2]]$Win.ID))))# 6736/10997=61.25%

unique(sort(c(intersect(Union.CANDf0.5_0.4_0.3[[2]]$Win.ID,Union.CANDf0.5_0.4_0.3[[3]]$Win.ID)))) #6736/10072=66.88%

unique(sort(c(intersect(Union.CANDf0.5_0.4_0.3[[6]]$Win.ID,Union.CANDf0.5_0.4_0.3[[7]]$Win.ID)))) #6368/9521=66.88%

unique(sort(c(intersect(Union.CANDf0.5_0.4_0.3[[7]]$Win.ID,Union.CANDf0.5_0.4_0.3[[6]]$Win.ID)))) #6368/9282=68.61%


#OUTLIER WINDOWS

unique(sort(c(intersect(Union.top0.5_0.4_0.3[[7]]$Win.ID,Union.top0.5_0.4_0.3[[2]]$Win.ID),  intersect(Union.top0.5_0.4_0.3[[7]]$Win.ID,Union.top0.5_0.4_0.3[[3]]$Win.ID), intersect(Union.top0.5_0.4_0.3[[7]]$Win.ID,Union.top0.5_0.4_0.3[[6]]$Win.ID)))) #1094/1270=86.14%


unique(sort(c(intersect(Union.top0.5_0.4_0.3[[3]]$Win.ID,Union.top0.5_0.4_0.3[[2]]$Win.ID),  intersect(Union.top0.5_0.4_0.3[[3]]$Win.ID,Union.top0.5_0.4_0.3[[6]]$Win.ID), intersect(Union.top0.5_0.4_0.3[[3]]$Win.ID,Union.top0.5_0.4_0.3[[7]]$Win.ID)))) #1057/1230=85.93%


unique(sort(c(intersect(Union.top0.5_0.4_0.3[[2]]$Win.ID,Union.top0.5_0.4_0.3[[3]]$Win.ID),  intersect(Union.top0.5_0.4_0.3[[2]]$Win.ID,Union.top0.5_0.4_0.3[[6]]$Win.ID), intersect(Union.top0.5_0.4_0.3[[2]]$Win.ID,Union.top0.5_0.4_0.3[[7]]$Win.ID)))) #1041/1248=83.41%


unique(sort(c(intersect(Union.top0.5_0.4_0.3[[6]]$Win.ID,Union.top0.5_0.4_0.3[[2]]$Win.ID),  intersect(Union.top0.5_0.4_0.3[[6]]$Win.ID,Union.top0.5_0.4_0.3[[3]]$Win.ID), intersect(Union.top0.5_0.4_0.3[[6]]$Win.ID,Union.top0.5_0.4_0.3[[7]]$Win.ID))))   #1093/1251=87.37%


#within continents

unique(sort(c(intersect(Union.top0.5_0.4_0.3[[3]]$Win.ID,Union.top0.5_0.4_0.3[[2]]$Win.ID))))# 924/1230=75.12%

unique(sort(c(intersect(Union.top0.5_0.4_0.3[[2]]$Win.ID,Union.top0.5_0.4_0.3[[3]]$Win.ID)))) #924/=1248=74.04%

unique(sort(c(intersect(Union.top0.5_0.4_0.3[[6]]$Win.ID,Union.top0.5_0.4_0.3[[7]]$Win.ID)))) #1008/1251=80.57%

unique(sort(c(intersect(Union.top0.5_0.4_0.3[[7]]$Win.ID,Union.top0.5_0.4_0.3[[6]]$Win.ID)))) #1008/1270=79.37%


#####

#for tf=0.5


length(unique(sort(c(intersect(CANDf0.5[[6]]$Win.ID,CANDf0.5[[2]]$Win.ID),  intersect(CANDf0.5[[6]]$Win.ID,CANDf0.5[[3]]$Win.ID), intersect(CANDf0.5[[6]]$Win.ID,CANDf0.5[[7]]$Win.ID)))))/nrow(CANDf0.5[[6]]) #75.12%
length(unique(sort(c(intersect(top829f0.5[[6]]$Win.ID,top829f0.5[[2]]$Win.ID),  intersect(top829f0.5[[6]]$Win.ID,top829f0.5[[3]]$Win.ID), intersect(top829f0.5[[6]]$Win.ID,top829f0.5[[7]]$Win.ID)))))/nrow(top829f0.5[[6]]) #85.64%


length(unique(sort(c(intersect(CANDf0.4[[6]]$Win.ID,CANDf0.4[[2]]$Win.ID),  intersect(CANDf0.4[[6]]$Win.ID,CANDf0.4[[3]]$Win.ID), intersect(CANDf0.4[[6]]$Win.ID,CANDf0.4[[7]]$Win.ID)))))/nrow(CANDf0.4[[6]]) #75.81%
length(unique(sort(c(intersect(top829f0.4[[6]]$Win.ID,top829f0.4[[2]]$Win.ID),  intersect(top829f0.4[[6]]$Win.ID,top829f0.4[[3]]$Win.ID), intersect(top829f0.4[[6]]$Win.ID,top829f0.4[[7]]$Win.ID)))))/nrow(top829f0.4[[6]]) #86.97


length(unique(sort(c(intersect(CANDf0.3[[6]]$Win.ID,CANDf0.3[[2]]$Win.ID),  intersect(CANDf0.3[[6]]$Win.ID,CANDf0.3[[3]]$Win.ID), intersect(CANDf0.3[[6]]$Win.ID,CANDf0.3[[7]]$Win.ID)))))/nrow(CANDf0.3[[6]]) #75.47%
length(unique(sort(c(intersect(top829f0.3[[6]]$Win.ID,top829f0.3[[2]]$Win.ID),  intersect(top829f0.3[[6]]$Win.ID,top829f0.3[[3]]$Win.ID), intersect(top829f0.3[[6]]$Win.ID,top829f0.3[[7]]$Win.ID)))))/nrow(top829f0.3[[6]]) #79.25



#Genes

#OUTLIER genes

unique(sort(c(intersect(TSI.top.genes[,1],LWK.top.genes[,1]),  intersect(TSI.top.genes[,1], YRI.top.genes[,1]), intersect(TSI.top.genes[,1],GBR.top.genes[,1])))) #203/231=88%

unique(sort(c(intersect(GBR.top.genes[,1],LWK.top.genes[,1]),  intersect(GBR.top.genes[,1], YRI.top.genes[,1]), intersect(GBR.top.genes[,1],TSI.top.genes[,1])))) #202/222=91%
unique(sort(c(intersect(LWK.top.genes[,1],YRI.top.genes[,1]),  intersect(LWK.top.genes[,1], GBR.top.genes[,1]), intersect(LWK.top.genes[,1],TSI.top.genes[,1])))) #208/249=83.5%
unique(sort(c(intersect(YRI.top.genes[,1],LWK.top.genes[,1]),  intersect(YRI.top.genes[,1], GBR.top.genes[,1]), intersect(YRI.top.genes[,1],TSI.top.genes[,1])))) #199/232=85.8%


#tf=0.5

unique(c(intersect(LWK.f0.5.genes[,1], YRI.f0.5.genes[,1]), intersect(LWK.f0.5.genes[,1], GBR.f0.5.genes[,1]), intersect(LWK.f0.5.genes[,1], TSI.f0.5.genes[,1]))) #917/1026
unique(c(intersect(LWK.f0.4.genes[,1], YRI.f0.4.genes[,1]), intersect(LWK.f0.4.genes[,1], GBR.f0.4.genes[,1]), intersect(LWK.f0.4.genes[,1], TSI.f0.4.genes[,1]))) #951/1099
unique(c(intersect(LWK.f0.3.genes[,1], YRI.f0.3.genes[,1]), intersect(LWK.f0.3.genes[,1], GBR.f0.3.genes[,1]), intersect(LWK.f0.3.genes[,1], TSI.f0.3.genes[,1]))) #1023/1182
#outlier

unique(c(intersect(LWK.f0.5.top829.genes[,1], YRI.f0.5.top829.genes[,1]), intersect(LWK.f0.5.top829.genes[,1], GBR.f0.5.top829.genes[,1]), intersect(LWK.f0.5.top829.genes[,1], TSI.f0.5.top829.genes[,1]))) #128/167


##
mean(c(length(intersect(LWK.top.genes[,1],YRI.top.genes[,1]))/nrow(LWK.top.genes),length(intersect(GBR.top.genes[,1],TSI.top.genes[,1]))/nrow(GBR.top.genes),length(intersect(YRI.top.genes[,1],LWK.top.genes[,1]))/nrow(YRI.top.genes),length(intersect(TSI.top.genes[,1],GBR.top.genes[,1]))/nrow(TSI.top.genes))) #78.8%

mean(c(length(intersect(LWK.cand.genes[,1],YRI.cand.genes[,1]))/nrow(LWK.cand.genes),length(intersect(GBR.cand.genes[,1],TSI.cand.genes[,1]))/nrow(GBR.cand.genes),length(intersect(YRI.cand.genes[,1],LWK.cand.genes[,1]))/nrow(YRI.cand.genes),length(intersect(TSI.cand.genes[,1],GBR.cand.genes[,1]))/nrow(TSI.cand.genes))) #79%



intersect(intersect(LWK.top.genes[,1], YRI.top.genes[,1]), intersect(GBR.top.genes[,1],TSI.top.genes[,1])) #102/249

mean(c(102/249, 102/232,102/222,102/231))
#[1] 0.4375779

#Cand

intersect(intersect(intersect(LWK.cand.genes[,1], YRI.cand.genes[,1]), GBR.cand.genes[,1]), TSI.cand.genes[,1]) #736

mean(rep(736,4)/c(nrow(LWK.cand.genes), nrow(YRI.cand.genes), nrow(GBR.cand.genes), nrow(TSI.cand.genes))) # 50%

