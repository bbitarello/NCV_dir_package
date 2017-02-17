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
Comparison=c(rep("G/(G+I)",8), rep("Ex/G",8), rep("Ex/all",8), rep("NSyn/Ex",8), rep("Reg/all",8),rep("Reg/(Reg+G)",8), rep("Reg/(Reg+I)",8), rep("Ex/(Ex+Reg)",8), rep("NSyn/(NSyn+Reg)",8)),
POP=rep(c("LWK", "LWK", "YRI","YRI","GBR","GBR", "TSI","TSI"),9), Set=rep(c("Significant","Outlier"),36)
)
#problem, we havep=0, so log(0)=-Inf. Need to do something about this. 
p<-c(rep(1,8),c(rep(0,6),0.003, 0),c(0.012,0, 0.013,0, 0.061,0.001, 0.088,0.001),c(0.002,0.009,0.004,0.022,0.037,0.115, 0.013,0.104),c(0.514,0.644,0.680,0.458,0.781,0.479,0.903,0.563),c(0,0.058,0,0.010,0.051,0.027,0.106,0.060),c(1,0.992,1,0.971,1, 0.950,1, 0.977),c(0.014, 0,0.005,0,0.016,0,0.012,0.003),c(0.771,0.026,0.563,0.036, 0.734,0.138,0.458,0.075))

#attention, need to do this so I don't get -Inf values!!
p[which(p==0)]<-0.001 
dt$logP<- abs(log(p))
factor(dt$POP, levels=c("LWK","YRI","GBR","TSI"))-> dt$POP
factor(dt$Set)-> dt$Set
factor(dt$Comparison, levels=c("G/(G+I)","Ex/G", "Ex/all","NSyn/Ex", "Reg/all", "Reg/(Reg+G)","Reg/(Reg+I)","Ex/(Ex+Reg)","NSyn/(NSyn+Reg)"))-> dt$Comparison

pdf('figures/test_enrich_plot_v1.pdf')
pdf('figures/test_enrich_plot_v2.pdf')
#this one works for sure
#ggplot(dt, aes(x=Comparison, y=logP, colour=POP, shape=Set)) + geom_point() + theme(axis.text.x=element_text(angle=45, hjust=1))

#works!!!
#ggplot(dt, aes(x=Comparison, y=logP, colour=POP, shape=Set)) + geom_point(size=3, alpha = 7/10) + scale_color_manual(values = wes_palette(n=4, name="Royal2")) + scale_shape_manual(values=c(17,19)) + theme(axis.text.x=element_text(angle=45, hjust=1))
c(abs(log(0.025)),abs(log(0.975)))->h

#v1
#ggplot(dt, aes(x=Comparison, y=logP, colour=POP, shape=Set)) + geom_point(size=2.5, alpha = 5/10, position = position_jitter(w = 0.5)) + scale_color_manual(values = c('RoyalBlue',"slateblue","darkolivegreen4", "darkolivegreen3")) + scale_shape_manual(values=c(17,19)) + theme(axis.text.x=element_text(angle=45, hjust=1)) + geom_hline(yintercept=c(abs(log(0.025)),abs(log(0.975))), color="darkgray", linetype="dashed") 

#v2
ggplot(dt, aes(x=Comparison, y=logP, colour=POP, shape=Set)) + geom_point(size=2.5, alpha = 5/10, position = position_jitter(w = 0.5)) + scale_color_manual(values = c('RoyalBlue',"slateblue","darkolivegreen4", "darkolivegreen3")) + scale_shape_manual(values=c(17,19)) + theme(axis.text.x=element_text(angle=45, hjust=1)) + geom_hline(yintercept=c(abs(log(0.025)),abs(log(0.975))), color="darkgray", linetype="dashed") +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))


dev.off()

