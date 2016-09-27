###############################################################################
#NCV scan for 1000g data
#dataAuthor: Barbara Bitarello
#BitarelloCreation: 17.12.2013
#Modified: 28.02.2014 by Cesare de Filippo
#Modified: 21.10.2014 by Barbara Bitarello
#BitarelloNote: attempting to fix two issues:run NCV for all pops at once and
#save window coordinates in output.
#Done
#Last Modified: 26.09.2016 by Barbara Bitarello 
###############################################################################
#NCV is calculated regardless of SNP density. Filtering per snp density should
#happen downstream in case we decide to change the threshold.

#in the future, modify this to have all pops run at once.
NCV.scan4<-function(INPUT.N, pop='YRI',FD=T, FD.N, WIN) {  
      ##INPUT.N : input data // ##FDs: fixed differences (human vs chimp
    ##reference)
        WIN[1,1]->beg; WIN[1,2]->end
        n<-100; if(pop=='PUR'){n<-88}  #it would be ideal to run for all pops at once and dothe rest of the commands on a per pop basis.
        if(FD==TRUE){ as.data.frame(FD.N)->z
        nifds<-dim(z)[1]}  #list of FDs between human and chimp
        nisnps<-dim(INPUT.N)[1] #number of SNPs in INPUT.N
      y2 <- as.data.frame(cbind(counts=as.numeric(INPUT.N[,pop])/n, pos=as.numeric(INPUT.N[,2]), ref=INPUT.N[,4], alt=INPUT.N[,5]), stringsAsFactors=F); #
       y2[,1]<-as.numeric(y2[,1]);y2[,2]<-as.numeric(y2[,2]);
	y3<-y2;
             #y3[,1] <- sapply(y2[,1], function(x) if (x>0.5){x<-1-x} else{x<-x})  #use minor allele frequency.
                #up until this point we have all original SNPs in y2. Now we
            #need to filter SNPs and FDs.
          if(FD==T){
	#real.snps <- y3[which(y3[,1] > 0),] 
        tmp<-which(z$pos %in% y2$pos) #positions present in SNP and FD file. Store them and then check what they are
                    tmp.vec<-NA;
                   if(length(tmp)){for (i in 1:length(tmp)){  #for each of these positions
    temp.pos<-z$pos[tmp[i]]; if(toupper(z$chimp[tmp[i]]) == toupper(y2[y2$pos==temp.pos,4])){#upper and lower case match.  
                                tmp.vec<-c(tmp.vec,tmp[i])}}  #store the positiuons which wuill be eliminated.
                            
    #if chimp in FD is != Alt in SNP, this is a SNP and also a FD. Keep both. 
    #no need to put any conditional statements for this
    tmp.vec<-tmp.vec[-1]  #eliminate the NA
	
    if(length(tmp.vec)){z2<-z[-tmp.vec,]}   #FD without the fake FDs.
        else{z2<-z} } else{z2<-z}  #check for positions present in SNPs that have f=1
                        #for alt allele and are absent from. FDs. Include them as FDs and exclude from SN                    
	 tmp.vec2<-0; tmp.vec3<-NA; tmp2<-which(y2$counts ==1);
            for (j in 1:length(tmp2)){  #for each of these positions
                temp.pos2<-y2$pos[tmp2[i]]; tmp.vec3<-c(tmp.vec3, temp.pos2); #exclude from SNPs
              if(!(temp.pos2 %in% z2$pos)){  #if this SNP position is not present in the FD
               tmp.vec2<-tmp.vec2+1}} # count number of FDs which should be included in NCV fd . ATTENTION: the FDinput bed will not be changed
	tmp.vec3<-tmp.vec3[-1]  #eliminate NA
                    fxdlen<-dim(z2)[1]+tmp.vec2 #if it is 0 no FDs will be added
                    if(length(tmp.vec3)){#if this vector has at least one position
                            real.snps2<-subset(y2, !(pos %in% tmp.vec3))}
                        if(length(tmp.vec3)==0){ #if there are no positions in tmp.vec3, all stays the same.
                                real.snps2<-y2}
                            real.snps3<-real.snps2[which(real.snps2$counts!=1),] #eliminate remaining snps with f=0 or f=1
                                polsites <- dim(real.snps3)[1] ;  #the 'real' number of SNPs used in NCV calculation.
				if(polsites>0){
                                tp<-as.numeric(c(real.snps3$counts,rep(0,fxdlen)));
				tp<-sapply(tp, function(x) if (x>0.5){x<-1-x} else{x<-x})} #use MAF
					else{tp<-c()}}
                ###############################################NCV without
                ###FD############################################
                    if(FD==F){
                            nifds<-0
                        real.snps4 <- y2[which(y2[,1] > 0 & y2[,1]<1),]  
                            polsites<-dim(real.snps4)[1];
tp<-as.numeric(real.snps4$counts); 
tp<-sapply(tp, function(x) if (x>0.5){x<-1-x} else{x<-x});
fxdlen<-0}
if(length(tp)>0){
ncdf1<-sqrt(sum((tp-0.1)^2)/(polsites+fxdlen)); 
ncdf2<-sqrt(sum((tp-0.2)^2)/(polsites+fxdlen));
ncdf3<-sqrt(sum((tp-0.3)^2)/(polsites+fxdlen));
ncdf4<-sqrt(sum((tp-0.4)^2)/(polsites+fxdlen)); ncdf5<-sqrt(sum((tp-0.5)^2)/(polsites+fxdlen));
                ################################################################fxdlen
          final<- data.frame(Beg.Win=beg, End.Win=end, Initial_seg_sites=nisnps, Initial_fds_sites=nifds, NCDf1=ncdf1,NCDf2=ncdf2, NCDf3=ncdf3,NCDf4=ncdf4, NCDf5=ncdf5,Nr.SNPs=polsites, Nr.FDs=fxdlen);
                return(final);}}

