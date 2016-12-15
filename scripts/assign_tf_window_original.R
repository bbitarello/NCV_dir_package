assign.tf<-function(df){
nrow(df)-> n
assigned.tf.per.window<-sapply(1:n, function(x) names(which.min(select(df[x,], Z.f0.5.P.val:Z.f0.3.P.val))))
min.p.value.per.window<-sapply(1:n, function(x) min(select(df[x,], Z.f0.5.P.val:Z.f0.3.P.val)))
cbind(assigned.ft=assigned.tf.per.window, min.p.value=min.p.value.per.window)-> res1
as.numeric(res1[,2])->res1[,2]
res1[which.min(res1[,2]),1]-> assigned.tf.gene
as.numeric(res1[which.min(res1[,2]),2])-> assigned.p.val.gene
return(list(assigned.per.window=res1,assigned.tf.per.gene= assigned.tf.gene, assigned.p.gene=assigned.p.val.gene))
}

