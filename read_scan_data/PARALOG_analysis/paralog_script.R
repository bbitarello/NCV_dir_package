###################################################################################
#	Bárbara Bitarello
#	Mast modified: 10.03.2016
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
#setDT(read.table('LWK_paralogs.txt', header=T, sep="\t"))-> LWK.cand.paralogs #reminder: this is not jsut LWK, but the shared sgnificant genes
#setDT(read.table('OR_LWK_paralogs2.txt', header=T, sep="\t"))-> OR.paralogs #reminder: this is the ORs just from the LWK Union CAND set...not sure why I did it like this.

nrow(ALL.GENES) #99818

length(unique(ALL.GENES$Associated_Gene_Name)) # 19,349

#read in list of scanned genes and filter:

fread('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/read_scan_data/bedfiles/all_scanned_genes.txt', header=F)-> scanned_genes

#filter

ALL.GENES  %>% dplyr::filter(Associated_Gene_Name %in% scanned_genes$V1) %>% as.data.table -> PAR_scanned
length(unique(PAR_scanned$Associated_Gene_Name)) #18593 so 40 genes don't show up that have weird annotation. That's ok.


fread('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/read_scan_data/bedfiles/prot.cod.Union.CAND.0.5_0.4_0.3_LWK.bed', header=F)$V -> LWK.signif


LWK.signif[grep("^OR",LWK.signif)][-27]-> ORs

PAR_scanned  %>% dplyr::filter(Associated_Gene_Name %in% LWK.signif$V1) %>% dplyr::filter(!(Associated_Gene_Name %in% ORs))

LWK.signif[grep("^OR",LWK.signif)][-27]-> ORs

PAR_scanned  %>% dplyr::filter(Associated_Gene_Name %in% LWK.signif) %>% dplyr::filter(!(Associated_Gene_Name %in% ORs))-> PAR_LWK_no_ORs

PAR_scanned %>% dplyr::filter(Associated_Gene_Name %in% ORs) -> PAR_ORs

na.omit(PAR_scanned)-> PAR_scanned.1
na.omit(PAR_LWK_no_ORs)-> PAR_LWK_no_ORs.1
na.omit(PAR_ORs) -> PAR_ORs.1


length(unique(PAR_scanned.1$Associated_Gene_Name)) #18593
length(unique(PAR_scanned.1$Associated_Gene_Name)) #12716 this mean 18593-12716=5877 genes have zero paralogues
length(unique(PAR_LWK_no_ORs$Associated_Gene_Name)) #1427
length(unique(PAR_LWK_no_ORs.1$Associated_Gene_Name)) # 1129
length(unique(PAR_ORs$Associated_Gene_Name)) #26
length(unique(PAR_ORs.1$Associated_Gene_Name)) #25

#factor(ALL.GENES[,1])-> ALL.GENES[,1] #13439 genes IDs

nrow(PAR_LWK_no_ORs) #8589

gcinfo(FALSE)
gc()


split(PAR_scanned.1, PAR_scanned.1$Associated_Gene_Name)-> big.list

split(PAR_LWK_no_ORs.1, PAR_LWK_no_ORs.1$Associated_Gene_Name)-> LWK.cand.paralogs.split
split(PAR_ORs.1, PAR_ORs.1$Associated_Gene_Name)-> OR.paralogs.split


#number os genes with zero paralogs at all:
#all scanned: 18593-12716=5877
#CAND: 1427-1129=298
#OR:26-25=1

PAR_scanned.1  %>% filter(!(Associated_Gene_Name %in% PAR_LWK_no_ORs.1$Associated_Gene_Name)) %>% group_by(Associated_Gene_Name) %>% filter(Chromosome.scaffold_name==Human_paralogue_chromosome.scaffold_name) %>% dplyr::summarise(N=n())-> Genomic

PAR_LWK_no_ORs.1 %>% group_by(Associated_Gene_Name) %>%  filter(Chromosome.scaffold_name==Human_paralogue_chromosome.scaffold_name) %>% dplyr::summarise(N=n())-> LWK.cand

PAR_ORs.1 %>% group_by(Associated_Gene_Name) %>% filter(Chromosome.scaffold_name==Human_paralogue_chromosome.scaffold_name) %>% dplyr::summarise(N=n())-> ORs


summary(Genomic$N) #min=1, max=76

X<-c(0,seq(1:76))

Gen<-vector('numeric',77)
#Signif<-vector('numeric',77)
Signif2<-vector('numeric',77)
OR_Genes<-vector('numeric',77)

c(0,as.numeric(labels(table(Genomic$N))[[1]]))-> tmp
for ( i in tmp[-1]){
y<-sum(table(Genomic$N)) #genes with paralogs in same chromosome #4,847 (exccuding the LWk candidates)
z<-sum(table((PAR_scanned.1  %>%  filter(!(Associated_Gene_Name %in% PAR_LWK_no_ORs.1$Associated_Gene_Name)) %>% group_by(Associated_Gene_Name) %>% summarise(N=n()))$N))-y  #genes with apralogs, but maybe not the same chr 11,587 (excl.LWK cand)
which(labels(table(Genomic$N))[[1]]==i)-> I
table(Genomic$N)[[I]]/(y+z)-> Gen[i+1]
}
z/(y+z)-> Gen[1] #almost 60% of scanned non candidate genes with at least one paralog anywhere have zero paralog in the same chromosome.
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

c(0,as.numeric(labels(table(LWK.cand$N))[[1]]))-> tmp4
for ( i in tmp4[-1]){
y<-sum(table(LWK.cand$N))
z<-sum(table((PAR_LWK_no_ORs.1 %>% group_by(Associated_Gene_Name) %>% summarise(N=n()))$N))-y
which(labels(table(LWK.cand$N))[[1]]==i)-> I
table(LWK.cand$N)[[I]]/(y+z)-> Signif2[i+1]
}
z/(z+y)-> Signif2[1]

c(0,as.numeric(labels(table(ORs$N))[[1]]))-> tmp3
for ( i in tmp3[-1]){
y<-sum(table(ORs$N))
z<-sum(table((PAR_ORs.1 %>% group_by(Associated_Gene_Name) %>% summarise(N=n()))$N))-y
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

