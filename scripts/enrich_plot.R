library(parallel)
library(SOAR)
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape)

devtools::install_github("karthik/wesanderson")
library(wesanderson)

dt<-data.table(
logP=rep(NA,72), 
Comparison=c(rep("Genic/(Genic+Intergenic)",8), rep("Exonic/Genic",8), rep("Exonic/all",8), rep("NSyn/Exonic",8), rep("Regulatory/all",8),rep("Regulatory/(Regulatory+Genic)",8), rep("Regulatory/(Regulatory+Intergenic)",8), rep("Exonic/(Exonic+Regulatory)",8), rep("NSyn/(NSyn+Regulatory)",8)),
POP=rep(c("LWK", "LWK", "YRI","YRI","GBR","GBR", "TSI","TSI"),9), Set=rep(c("Significant","Outlier"),36)
)
#problem, we havep=0, so log(0)=-Inf. Need to do something about this. 
p<-c(rep(1,8),c(rep(0,6),0.003, 0),c(0.012,0, 0.013,0, 0.061,0.001, 0.088,0.001),c(0.002,0.009,0.004,0.022,0.037,0.115, 0.013,0.104),c(0.514,0.644,0.680,0.458,0.781,0.479,0.903,0.563),c(0,0.058,0,0.010,0.051,0.027,0.106,0.060),c(1,0.992,1,0.971,1, 0.950,1, 0.977),c(0.014, 0,0.005,0,0.016,0,0.012,0.003),c(0.771,0.026,0.563,0.036, 0.734,0.138,0.458,0.075))

#attention, need to do this so I don't get -Inf values!!
p[which(p==0)]<-0.001 
dt$logP<- abs(log(p))
factor(dt$POP, levels=c("LWK","YRI","GBR","TSI"))-> dt$POP
factor(dt$Set)-> dt$Set
factor(dt$Comparison, levels=c("Exonic/all","Exonic/Genic", "NSyn/Exonic", "Regulatory/all", "Exonic/(Exonic+Regulatory)","NSyn/(NSyn+Regulatory)","Regulatory/(Regulatory+Genic)","Regulatory/(Regulatory+Intergenic)","Genic/(Genic+Intergenic)"))-> dt$Comparison
dt %>% dplyr::filter(Comparison %in% c("Exonic/all","Exonic/Genic", "NSyn/Exonic", "Regulatory/all", "Exonic/(Exonic+Regulatory)","NSyn/(NSyn+Regulatory)"))-> dt2
#pdf('figures/test_enrich_plot_v1.pdf')
#png('figures/test_enrich_plot_v2.png')

tiff(file="figures/enrich2.tiff", height = 90.5, width = 150, units = 'mm', compression = "lzw", res = 500)
h<-log(0.05);a<-expression(paste0('-log(', italic('p'), ')'))
#ggplot(dt2, aes(x=Comparison, y=logP, colour=POP, shape=Set)) + geom_point(size=2.5, alpha = 5/10, position = position_jitter(w = 0.5)) + scale_color_manual(values = c('RoyalBlue',"slateblue","darkolivegreen4", "darkolivegreen3")) + scale_shape_manual(values=c(17,19)) + theme(axis.text.x=element_text(angle=45, hjust=1), axis.title.x=element_blank()) + geom_hline(yintercept=abs(log(0.05)), color="darkgray", linetype="dashed") + ylab(expression(-log(italic(p))))
#v2
ggplot(dt2, aes(x=Comparison, y=logP, colour=POP, shape=Set)) + geom_point(size=2.5, alpha = 5/10, position = position_jitter(w = 0.5)) + scale_color_manual(values = c('RoyalBlue',"slateblue","darkolivegreen4", "darkolivegreen3")) + scale_shape_manual(values=c(17,19)) + theme_bw() +theme(axis.text.x=element_text(angle=45, hjust=1), plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme(panel.border= element_blank()) + theme(axis.line.y = element_line(color="black", size = 0.5),axis.title.x=element_blank(), axis.line.x = element_line(color="black", size = 0.5)) + geom_hline(yintercept=abs(log(0.05)), color="darkgray", linetype="dashed") + ylab(expression(-log(italic(p)))) +  coord_cartesian(ylim=c(00,8))
dev.off()

