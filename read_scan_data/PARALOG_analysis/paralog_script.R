###################################################################################
#	Bárbara Bitarello
#	Mast modified: 24.11.2016
#	This  is based on the old scripts, but hopefully with some improvements.
#
####################################################################################

library(parallel)
library(SOAR)
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape)

setDT(read.table('Background.txt', header=T, sep="\t"))-> ALL.GENES
setDT(read.table('LWK_paralogs.txt', header=T, sep="\t"))-> LWK.cand.paralogs
setDT(read.table('OR_LWK_paralogs2.txt', header=T, sep="\t"))-> OR.paralogs

nrow(ALL.GENES) #69108

length(unique(ALL.GENES$Associated_Gene_Name)) # 13447
na.omit(ALL.GENES)-> ALL.GENES.1
#factor(ALL.GENES[,1])-> ALL.GENES[,1] #13439 genes IDs

length(unique(ALL.GENES.1$Associated_Gene_Name)) #9224  #these are the genes left after na.omit
nrow(ALL.GENES.1) #64871
remove(ALL.GENES)
gc()


nrow(LWK.cand.paralogs) #8832
na.omit(LWK.cand.paralogs)-> LWK.cand.paralogs1
remove(LWK.cand.paralogs); gc()
nrow(LWK.cand.paralogs1) #8528
nrow(OR.paralogs) #235
na.omit(OR.paralogs)-> OR.paralogs1
remove(OR.paralogs); gc()
nrow(OR.paralogs1) #234

gcinfo(FALSE)
gc()


split(ALL.GENES.1, ALL.GENES.1$Associated_Gene_Name)-> big.list
split(LWK.cand.paralogs1, LWK.cand.paralogs1$Associated_Gene_Name-> LWK.cand.paralogs.split
split(OR.paralogs1, OR.paralogs1$Ensembl_Associated_Gene_Name)-> OR.paralogs.split



ALL.GENES.1 %>% group_by(Associated_Gene_Name) %>% filter(Chromosome_Name==Human_Paralogue_Chromosome_Name) %>% summarise(N=n())-> Genomic
LWK.cand.paralogs1 %>% group_by(Associated_Gene_Name) %>%  filter(Chromosome_Name==Human_Paralogue_Chromosome_Name) %>% summarise(N=n())-> LWK.cand
LWK.cand.paralogs1 %>% filter(!(Associated_Gene_Name %in% OR.paralogs1$Associated_Gene_Name)) %>% group_by(Associated_Gene_Name) %>% filter(Chromosome_Name==Human_Paralogue_Chromosome_Name) %>%summarise(N=n()) -> LWK.no.OR
OR.paralogs1 %>% group_by(Associated_Gene_Name) %>% summarise(N=n())-> ORs


summary(Genomic$N) #min=1, max=80

X<-seq(1:80)

Gen<-vector('numeric',80)
Signif<-vector('numeric',80)
Signif2<-vector('numeric',80)
OR_Genes<-vector('numeric',80)

as.numeric(labels(table(Genomic$N))[[1]])-> tmp
for ( i in tmp){
y<-sum(table(Genomic$N))
which(labels(table(Genomic$N))[[1]]==i)-> I
table(Genomic$N)[[I]]/y-> Gen[i]
}


as.numeric(labels(table(LWK.cand$N))[[1]])-> tmp2
for ( i in tmp2){
y<-sum(table(LWK.cand$N))
which(labels(table(LWK.cand$N))[[1]]==i)-> I
table(LWK.cand$N)[[I]]/y-> Signif[i]
}


as.numeric(labels(table(LWK.no.OR$N))[[1]])-> tmp4
for ( i in tmp4){
y<-sum(table(LWK.no.OR$N))
which(labels(table(LWK.no.OR$N))[[1]]==i)-> I
table(LWK.no.OR$N)[[I]]/y-> Signif2[i]
}


as.numeric(labels(table(ORs$N))[[1]])-> tmp3
for ( i in tmp3){
y<-sum(table(ORs$N))
which(labels(table(ORs$N))[[1]]==i)-> I
table(ORs$N)[[I]]/y-> OR_Genes[i]
}


to_plot <- data.frame(x=X,Genomic=Gen,Significant=Signif,Significant_no_ORs=Signif2, ORs=OR_Genes)
to_plot<-to_plot[1:40,] #after 38 there is only 3 genes in the background list
melted<-melt(to_plot, id="x")
as.factor(melted$x)-> melted$x




pdf('paralogs_LWK_same_chr.pdf')
print(ggplot(melted,aes(x=x,y=value,fill=variable)) + geom_bar(stat="identity", alpha=.4) + scale_fill_manual(values=c('darkgray', 'darkolivegreen', 'steelblue', 'darkmagenta')) + ylab("Proportion of genes") + xlab("# of paralogs on the same chromosome")+theme(panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), axis.text.x=element_text(size=10, angle=45)) + scale_x_discrete(breaks=seq(from=1, to=50, by=5)))

dev.off()

