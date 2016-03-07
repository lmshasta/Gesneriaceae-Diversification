#geosse script for Andes vs. rest of globe

library(diversitree)

pruned.trees<-read.nexus("pruned.trees.tre")

data=read.table("data.africa.txt", row.names=1)
africa<-data[,1]
names(africa)<-row.names(data)

#skeletal tree approach for random sampling
sampling.f<-c(1.0,0.75,0.25)

#full bisse model
fit<-lapply(1:length(pruned.trees), function(x){
  model<-make.geosse(tree=pruned.trees[[x]], states=africa, sampling.f=sampling.f)
  start.values<-starting.point.geosse(pruned.trees[[x]])
  res<-find.mle(model, x.init=start.values, method="subplex")
  cat("Job #", x, "finished\n", sep="")
  return(res)
})

fit.est.param.table<-as.data.frame(do.call("rbind", 
                                           lapply(fit, function(x)( c(x$par, npar=length(x$par), lnl=x$lnLik)))))

#Write the table to a file
output<-paste(Sys.getenv("HOME"), "/fit.txt", sep='')
write.table(fit.est.param.table, file=output, row.names=F, quote=F, sep="\t")

#constrain same speciation within-regions
fit.1<-lapply(1:length(pruned.trees), function(x){
  model<-make.geosse(tree=pruned.trees[[x]], states=africa, sampling.f=sampling.f)
  lik.1<-constrain(model, sA ~ sB, sAB ~ 0) 
  start.values<-starting.point.geosse(pruned.trees[[x]])
  res.1<-find.mle(lik.1, start.values[argnames(lik.1)])
  cat("job #", x, "finished\n", sep="")
  return(res.1)
})

fit.1.est.param.table<-as.data.frame(do.call("rbind", 
                                             lapply(fit.1, function(x)( c(x$par, npar=length(x$par), lnl=x$lnLik)))))

#write the table to a file
output<-paste(Sys.getenv("HOME"), "/fit.1.txt", sep='')
write.table(fit.1.est.param.table, file=output, row.names=F, quote=F, sep="\t")

#constrain between-region speciation
fit.2<-lapply(1:length(pruned.trees), function(x){
  model<-make.geosse(tree=pruned.trees[[x]], states=africa, sampling.f=sampling.f)
  lik.2<-constrain(model, sAB ~ 0)
  start.values<-starting.point.geosse(pruned.trees[[x]])
  res.2<-find.mle(lik.2, start.values[argnames(lik.2)])
  cat("job #", x, "finished\n", sep="")
  return(res.2)
})

fit.2.est.param.table<-as.data.frame(do.call("rbind", 
                                             lapply(fit.2, function(x)( c(x$par, npar=length(x$par), lnl=x$lnLik)))))

#write the table to a file
output<-paste(Sys.getenv("HOME"), "/fit.2.txt", sep='')
write.table(fit.2.est.param.table, file=output, row.names=F, quote=F, sep="\t")

#constrain dispersal between regions
fit.3<-lapply(1:length(pruned.trees), function(x){
  model<-make.geosse(tree=pruned.trees[[x]], states=africa, sampling.f=sampling.f)
  lik.3<-constrain(model, dA ~ dB)
  start.values<-starting.point.geosse(pruned.trees[[x]])
  res.3<-find.mle(lik.3, start.values[argnames(lik.3)])
  cat("job #", x, "finished\n", sep="")
  return(res.3)
})

fit.3.est.param.table<-as.data.frame(do.call("rbind", 
                                             lapply(fit.3, function(x)( c(x$par, npar=length(x$par), lnl=x$lnLik)))))

#write the table to a file
output<-paste(Sys.getenv("HOME"), "/fit.3.txt", sep='')
write.table(fit.3.est.param.table, file=output, row.names=F, quote=F, sep="\t")

#constrain within-region extinction
fit.4<-lapply(1:length(pruned.trees), function(x){
  model<-make.geosse(tree=pruned.trees[[x]], states=africa, sampling.f=sampling.f)
  lik.4<-constrain(model, xA ~ xB)
  start.values<-starting.point.geosse(pruned.trees[[x]])
  res.4<-find.mle(lik.4, start.values[argnames(lik.4)])
  cat("job #", x, "finished\n", sep="")
  return(res.4)
})

fit.4.est.param.table<-as.data.frame(do.call("rbind", 
                                             lapply(fit.4, function(x)( c(x$par, npar=length(x$par), lnl=x$lnLik)))))

#write the table to a file
output<-paste(Sys.getenv("HOME"), "/fit.4.txt", sep='')
write.table(fit.4.est.param.table, file=output, row.names=F, quote=F, sep="\t")
