#######
#	10.12.2015
#	Last modified: 15.11.2016
#
#######


#proportions of scanned windows per chromosome, vs simulation based vs top 816 (YRI)


library(ggplot2)
library(reshape)



#LUHYA (LWK)

#genomic

chrs.gen<-unlist(mclapply(1:22, function(x) nrow(filter(list.SCAN[[2]], Chr==x))))
x<-seq(1:22)

chrs.cand<-unlist(mclapply(1:22, function(x) nrow(filter(Union.CANDf0.5_0.4_0.3[[2]], Chr==paste0('chr',x)))))

chrs<-unlist(mclapply(1:22, function(x) nrow(filter(Union.top0.5_0.4_0.3[[2]], Chr==paste0('chr',x)))))

y1<-chrs.gen/sum(chrs.gen)

y2<-chrs.cand/sum(chrs.cand)

y3<-chrs/sum(chrs)



to_plot <- data.frame(x=x,Genomic=y1,Significant=y2, Outliers=y3)
melted<-melt(to_plot, id="x")
as.factor(melted$x)-> melted$x

pdf('figures/chromosome.distrib.windows.LWK.pdf')
print(ggplot(melted,aes(x=x,y=value,fill=variable)) + geom_bar(stat="identity", alpha=.3) + scale_fill_manual(values=c('darkgray', 'darkolivegreen', 'steelblue')) + ylab("Proportion of windows") + xlab("Chromosome")+ theme(panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank()))
dev.off()

