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
setDT(read.table('LWK_paralogs.txt', header=T, sep="\t"))-> LWK.cand.paralogs #reminder: this is not jsut LWK, but the shared sgnificant genes
setDT(read.table('OR_LWK_paralogs2.txt', header=T, sep="\t"))-> OR.paralogs #reminder: this is the ORs just from the LWK Union CAND set...not sure why I did it like this.

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
setDT(LWK.cand.paralogs1 %>% filter(!(Associated_Gene_Name %in% OR.paralogs$Associated_Gene_Name)))-> no.ORs

nrow(OR.paralogs) #283
na.omit(OR.paralogs)-> OR.paralogs1
remove(OR.paralogs); gc()
nrow(OR.paralogs1) #282

gcinfo(FALSE)
gc()


split(ALL.GENES.1, ALL.GENES.1$Associated_Gene_Name)-> big.list
split(LWK.cand.paralogs1, LWK.cand.paralogs1$Associated_Gene_Name)-> LWK.cand.paralogs.split
split(OR.paralogs1, OR.paralogs1$Associated_Gene_Name)-> OR.paralogs.split


#number os genes with zero paralogs:
#all scanned: 4,223
#CAND 302
#OR:1
ALL.GENES.1 %>% filter(!(Associated_Gene_Name %in% LWK.cand.paralogs1$Associated_Gene_Name)) %>% group_by(Associated_Gene_Name) %>% filter(Chromosome_Name==Human_Paralogue_Chromosome_Name) %>% summarise(N=n())-> Genomic
LWK.cand.paralogs1 %>% group_by(Associated_Gene_Name) %>%  filter(Chromosome_Name==Human_Paralogue_Chromosome_Name) %>% summarise(N=n())-> LWK.cand
LWK.cand.paralogs1 %>% filter(!(Associated_Gene_Name %in% OR.paralogs1$Associated_Gene_Name)) %>% group_by(Associated_Gene_Name) %>% filter(Chromosome_Name==Human_Paralogue_Chromosome_Name) %>%summarise(N=n()) -> LWK.no.OR
OR.paralogs1 %>% group_by(Associated_Gene_Name) %>% summarise(N=n())-> ORs


summary(Genomic$N) #min=1, max=76

X<-c(0,seq(1:76))

Gen<-vector('numeric',77)
#Signif<-vector('numeric',77)
Signif2<-vector('numeric',77)
OR_Genes<-vector('numeric',77)

c(0,as.numeric(labels(table(Genomic$N))[[1]]))-> tmp
for ( i in tmp[-1]){
y<-sum(table(Genomic$N))
z<-sum(table((ALL.GENES.1 %>% group_by(Associated_Gene_Name) %>% summarise(N=n()))$N))-y
which(labels(table(Genomic$N))[[1]]==i)-> I
table(Genomic$N)[[I]]/(y+z)-> Gen[i+1]
}
z/(y+z)-> Gen[1]
remove(y,z)

#c(0,as.numeric(labels(table(LWK.cand$N))[[1]]))-> tmp2
#for ( i in tmp2[-1]){
#y<-sum(table(LWK.cand$N))
#z<-sum(table((LWK.cand.paralogs1 %>% group_by(Associated_Gene_Name) %>% summarise(N=n()))$N))-y
#which(labels(table(LWK.cand$N))[[1]]==i)-> I
#table(LWK.cand$N)[[I]]/(y+z)-> Signif[i+1]
#}
#z/(y+z)->Signif[1]
#remove(y,z)

c(0,as.numeric(labels(table(LWK.no.OR$N))[[1]]))-> tmp4
for ( i in tmp4[-1]){
y<-sum(table(LWK.no.OR$N))
z<-sum(table((no.ORs %>% group_by(Associated_Gene_Name) %>% summarise(N=n()))$N))-y
which(labels(table(LWK.no.OR$N))[[1]]==i)-> I
table(LWK.no.OR$N)[[I]]/(y+z)-> Signif2[i+1]
}
z/(z+y)-> Signif2[1]

c(0,as.numeric(labels(table(ORs$N))[[1]]))-> tmp3
for ( i in tmp3[-1]){
y<-sum(table(ORs$N))
z<-sum(table((OR.paralogs1 %>% group_by(Associated_Gene_Name) %>% summarise(N=n()))$N))-y
which(labels(table(ORs$N))[[1]]==i)-> I
table(ORs$N)[[I]]/(y+z)-> OR_Genes[i+1]
}
z/(y+z)-> OR_Genes[1]
remove(y,z)

#to_plot <- data.frame(x=X,Genomic=Gen,Significant=Signif,Significant_no_ORs=Signif2, ORs=OR_Genes)

to_plot <- data.frame(x=X,Genomic=Gen,Significant_no_ORs=Signif2, ORs=OR_Genes)

to_plot<-to_plot[0:40,] #after 38 there is only 3 genes in the background list
melted<-melt(to_plot, id="x")
as.factor(melted$x)-> melted$x


pdf('for_paper_paralogs.pdf')

print(ggplot(melted, aes(x=x, y=value))+ geom_bar(aes(fill = variable), position = "dodge", stat="identity", alpha=.4) + scale_fill_manual(values=c('darkgray', 'darkolivegreen', 'darkmagenta'))  + ylab("Proportion of genes") + xlab("# of paralogs on the same chromosome") + theme(panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), axis.text.x=element_text(size=10, angle=45)) + scale_x_discrete(breaks=seq(from=0, to=79, by=5)))

dev.off()

#pdf('paralogs_LWK_same_chr.pdf')

pdf('paralogs_LWK_same_chr1.pdf')

print(ggplot(melted,aes(x=x,y=value,fill=variable)) + geom_bar(stat="identity", alpha=.4) + scale_fill_manual(values=c('darkgray', 'darkolivegreen', 'steelblue', 'darkmagenta')) + ylab("Proportion of genes") + xlab("# of paralogs on the same chromosome")+theme(panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), axis.text.x=element_text(size=10, angle=45)) + scale_x_discrete(breaks=seq(from=0, to=49, by=5)))

dev.off()

