############################################################################################################################
#
#       Barbara D Bitarello
#	Created:02.10.2016
#       Last modified: 02.10.2016
#	A script to start processing scan data.
#
#################################################################################################################################


#first, load the scan data


#load packages

library(parallel)
library(SOAR)
library(ggplot2)
Sys.setenv(R_LOCAL_CACHE="estsession")

load('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/read_scan_data/All.Results.Final.RData')

##########################skip the next block, as it has already been saved in the R object ############################
pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")

names(All.Results.Final)<-pops

#obsolete because we don't use the cov info anymore.#

#cov.win<-read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/windows_coordinates_cov.bed.gz')

#names(cov.win)<-c('CHR', 'Beg.Win', 'End.Win', 'Nr.Map.Seg', 'Total.Cov.Leng', 'Total.Win.Leng', 'Proportion.Covered')

#add coverage to dataframes

#cov.win[order(cov.win$CHR, cov.win$Beg.Win),]-> cov.win2   #the coverage values are not in oder (the windows)

#now windows are sorted in NCV output. 

#remove(cov.win)

#mclapply(All.Results.Final, function(x) cbind(x, Proportion.Covered=cov.win2$Proportion.Covered))->All.Results.Final

#objectName<-'All.Results.Final'

#save(list=objectName, file= 'All.Results.Final.RData')   #this workspace is already saved in the directory.

#remove(All.Results.Final)
                                                                                                                        

#load('/mnt/sequencedb/PopGen/barbara/scan_may_2014/All.Results.Final.RData')


#starts here

mclapply(All.Results.Final, function(x) cbind(x,Nr.IS=(x$Nr.SNPs+x$Nr.FDs)))-> All.Res

#now we have a data frame which contain the Nr.US

mclapply(All.Res, function(x) cbind(x,PtoD=x$Nr.SNPs/(x$Nr.FDs+1)))-> All.Res1

mclapply(All.Res1, function(x) dim(x)) #dimensions of initial data

# 1705970 #update:1,543,026 for YRI. But actually it varies between pops now. WHY???
#update 05.10/2016: between 1532463 (FIN) and 1543712 (AWS)
#update 07.10.2016:1551642 for all pops!
#mclapply(All.Res1, function(x) subset(x, Nr.IS>=4))->All.Res.4.IS
#yri:1698893 #update:


#New: add WIn.ID column

mclapply(1:7, function(x) with(All.Res1[[x]], paste0(Chr, "|", Beg.Win, "|", End.Win)))-> tmp

mclapply(1:7, function(x) cbind(All.Res1[[x]], Win.ID=tmp[[x]],stringsAsFactors=FALSE))-> tmp2

tmp2->All.Res1
remove(tmp2)


#exploring number of IS (new 07/10/2016)


#YRI

Nr.IS<-as.factor(with(All.Results.Final[[3]], Nr.SNPs+Nr.FDs))
temp<-data.frame(NCD2=All.Results.Final[[3]]$NCDf5, IS=Nr.IS)
temp.df<-melt(temp, id.vars='IS')

#99.4 % of all windows have <=100 informative sites so i don't need to plot everything.
temp.df3<-subset(temp.df, as.numeric(IS)>=19 & as.numeric(IS)<=100)
temp.df2<-subset(temp.df, as.numeric(IS)<=100)
pdf('violin.Nr.IS.YRI.empirical.pdf')
ggplot(temp.df2, aes(IS, value))+scale_x_discrete(breaks = seq(1,100, by=5)) +geom_violin() + stat_summary(fun.y=median, geom="point", color='cornflowerblue', size=1) +  ylab("NCD2") + xlab("Informative Sites") + geom_vline(xintercept = 10, colour="orange", linetype = "longdash")

dev.off()


#LWK


Nr.IS<-as.factor(with(All.Results.Final[[2]], Nr.SNPs+Nr.FDs))
temp<-data.frame(NCD2=All.Results.Final[[2]]$NCDf5, IS=Nr.IS)
temp.df<-melt(temp, id.vars='IS')

#99.3 % of all windows have <=100 informative sites so i don't need to plot everything.
temp.df3<-subset(temp.df, as.numeric(IS)>=19 & as.numeric(IS)<=100)
temp.df2<-subset(temp.df, as.numeric(IS)<=100)
pdf('violin.Nr.IS.LWK.empirical.pdf')
ggplot(temp.df2, aes(IS, value))+scale_x_discrete(breaks = seq(1,100, by=5)) +geom_violin() + stat_summary(fun.y=median, geom="point", color='cornflowerblue', size=1) +  ylab("NCD2") + xlab("Informative Sites")  + geom_vline(xintercept = 10, colour="orange", linetype = "longdash")
dev.off()


#GBR

Nr.IS<-as.factor(with(All.Results.Final[[6]], Nr.SNPs+Nr.FDs))
temp<-data.frame(NCD2=All.Results.Final[[6]]$NCDf5, IS=Nr.IS)
temp.df<-melt(temp, id.vars='IS')

#99.6 % of all windows have <=100 informative sites so i don't need to plot everything.
temp.df3<-subset(temp.df, as.numeric(IS)>=15 & as.numeric(IS)<=100)
temp.df2<-subset(temp.df, as.numeric(IS)<=100)

pdf('violin.Nr.IS.GBR.empirical.pdf')

ggplot(temp.df2, aes(IS, value))+scale_x_discrete(breaks = seq(1,100, by=5)) +geom_violin() + stat_summary(fun.y=median, geom="point", color='cornflowerblue', size=1) +  ylab("NCD2") + xlab("Informative Sites")  + geom_vline(xintercept = 10, colour="orange", linetype = "longdash")

dev.off()


#TSI

Nr.IS<-as.factor(with(All.Results.Final[[7]], Nr.SNPs+Nr.FDs))
temp<-data.frame(NCD2=All.Results.Final[[7]]$NCDf5, IS=Nr.IS)
temp.df<-melt(temp, id.vars='IS')

#99.7% of all windows have <=100 informative sites so i don't need to plot everything.
temp.df3<-subset(temp.df, as.numeric(IS)>=19 & as.numeric(IS)<=100)
temp.df2<-subset(temp.df, as.numeric(IS)<=100)
pdf('violin.Nr.IS.TSI.empirical.pdf')
ggplot(temp.df2, aes(IS, value))+scale_x_discrete(breaks = seq(1,100, by=5)) +geom_violin() + stat_summary(fun.y=median, geom="point", color='cornflowerblue', size=1) +  ylab("NCD2") + xlab("Informative Sites")  + geom_vline(xintercept = 10, colour="orange", linetype = "longdash")

dev.off()




#will need to check if this filter is still appropriate:
#update: i think >=9 is appropriate for all pops now (see violin plots)

#this filter would result in about 20,000 windows exluded per pop.

unique(sort(
c(All.Res1[[3]]$Win.ID[which(All.Res1[[3]]$Nr.IS<9)], 
All.Res1[[2]]$Win.ID[which(All.Res1[[2]]$Nr.IS<9)], 
All.Res1[[1]]$Win.ID[which(All.Res1[[1]]$Nr.IS<9)], 
All.Res1[[4]]$Win.ID[which(All.Res1[[4]]$Nr.IS<9)], 
All.Res1[[5]]$Win.ID[which(All.Res1[[5]]$Nr.IS<9)],
All.Res1[[6]]$Win.ID[which(All.Res1[[6]]$Nr.IS<9)], 
All.Res1[[7]]$Win.ID[which(All.Res1[[7]]$Nr.IS<9)]
)))->filt.Nr.IS #or 26228 in total


mclapply(All.Res1, function(x) subset(x, !(Win.ID %in% filt.Nr.IS)))->All.Res2

#which results in

unlist(lapply(All.Res2, function(x) nrow(x)))  #1,525,424 windows per pop, i.e, only 1.7% os fara was removed with the IS filter. Totally worth it.
#mclapply(All.Res2, function(x) subset(x, Proportion.Covered>=0.5))->All.Res.filtered

#lapply(All.Res.4.IS, function(x) cbind(x, P.val.NCVf0.5=rep(NA, dim(x)[1]), P.val.NCVf0.4=rep(NA, dim(x)[1]), P.val.NCVf0.3=rep(NA, dim(x)[1]), P.val.NCVf0.2=rep(NA, dim(x)[1]), P.val.NCVf0.1=rep(NA, dim(x)[1])))->All.Res.4.IS

#lapply(All.Res.4.IS, function(x) subset(x, Proportion.Covered>=0.5))-> All.Res.4.IS.prop50
#dim YRI:1627870 

#objectName<-'All.Res.4.IS.prop50'
#objectName2<-'All.Res.filtered'
setwd('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/read_scan_data/')

#save(list=objectName, file= 'All.Res.4.IS.prop50.RData')
#save(list=objectName2, file='All.Res.filtered.RData')

#Store(All.Res.4.IS.prop50)
#Store(All.Res.filtered)

Store(All.Res2)
Objects()

objectName<-'All.Res2'

save(list=objectName, file= 'Results.After.IS.filter.RData')
###################
###################
###################