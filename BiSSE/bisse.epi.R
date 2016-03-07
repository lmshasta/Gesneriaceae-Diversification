library(diversitree)


#import phylogenetic tress
trees<-read.tree("treefile")

#import dataframe
data=read.table("datafile", row.names=1)
epi<-data[,1]
names(epi)<-row.names(data)

#skeletal tree approach for random sampling
sampling.f<-c(0.267,0.198)

#full bisse model
fit<-lapply(1:length(pruned.trees), function(x){
  model<-make.bisse(tree=pruned.trees[[x]], states=epi, sampling.f=sampling.f)
  start.values<-starting.point.bisse(pruned.trees[[x]],yule=TRUE)
  res<-find.mle(model, x.init=start.values, method="subplex")
  cat("Job #", x, "finished\n", sep="")
  return(res)
})

fit.est.param.table<-as.data.frame(do.call("rbind", 
      lapply(fit, function(x)( c(x$par, npar=length(x$par), lnl=x$lnLik)))))

#Write the table to a file
output<-paste(Sys.getenv("HOME"), "/fit.txt", sep='')
write.table(fit.est.param.table, file=output, row.names=F, quote=F, sep="\t")

#constrain lambda
fit.1<-lapply(1:length(pruned.trees), function(x){
  model<-make.bisse(tree=pruned.trees[[x]], states=epi, sampling.f=sampling.f)
  lik.1<-constrain(model, lambda1 ~ lambda0) 
  start.values<-starting.point.bisse(pruned.trees[[x]])
  res.1<-find.mle(lik.1, start.values[argnames(lik.1)])
  cat("job #", x, "finished\n", sep="")
  return(res.1)
})

fit.1.est.param.table<-as.data.frame(do.call("rbind", 
      lapply(fit.1, function(x)( c(x$par, npar=length(x$par), lnl=x$lnLik)))))

#write the table to a file
output<-paste(Sys.getenv("HOME"), "/fit.1.txt", sep='')
write.table(fit.1.est.param.table, file=output, row.names=F, quote=F, sep="\t")

#mcmc - equal extinction
  model<-make.bisse(tree=pruned.trees[[1]],states=epi,sampling.f=sampling.f)
  lik<-constrain(model,mu1 ~ mu0)
  start.values<-starting.point.bisse(pruned.trees[[1]])
  res<-find.mle(model,x.init=start.values[argnames(lik)],method="subplex",control=list(maxit=50000))
  set.seed(1)
  prior<-make.prior.exponential(1/(2*(res$par[2]-res$par[3])))
  samples<-mcmc(lik,res$par[argnames(lik)],nsteps=1000000,w=.3,prior=prior,print.every=100)
}

#constrain mu
fit.2<-lapply(1:length(pruned.trees), function(x){
  model<-make.bisse(tree=pruned.trees[[x]], states=epi, sampling.f=sampling.f)
  lik.2<-constrain(model, mu1 ~ mu0)
  start.values<-starting.point.bisse(pruned.trees[[x]])
  res.2<-find.mle(lik.2, start.values[argnames(lik.2)])
  cat("job #", x, "finished\n", sep="")
  return(res.2)
})

fit.2.est.param.table<-as.data.frame(do.call("rbind", 
        lapply(fit.2, function(x)( c(x$par, npar=length(x$par), lnl=x$lnLik)))))

#write the table to a file
output<-paste(Sys.getenv("HOME"), "/fit.2.txt", sep='')
write.table(fit.2.est.param.table, file=output, row.names=F, quote=F, sep="\t")

#constrain q
fit.3<-lapply(1:length(pruned.trees), function(x){
  model<-make.bisse(tree=pruned.trees[[x]], states=epi, sampling.f=sampling.f)
  lik.3<-constrain(model, q01 ~ q10)
  start.values<-starting.point.bisse(pruned.trees[[x]])
  res.3<-find.mle(lik.3, start.values[argnames(lik.3)])
  cat("job #", x, "finished\n", sep="")
  return(res.3)
})

fit.3.est.param.table<-as.data.frame(do.call("rbind", 
     lapply(fit.3, function(x)( c(x$par, npar=length(x$par), lnl=x$lnLik)))))

#write the table to a file
output<-paste(Sys.getenv("HOME"), "/fit.3.txt", sep='')
write.table(fit.3.est.param.table, file=output, row.names=F, quote=F, sep="\t")

#no reversal
fit.4<-lapply(1:length(pruned.trees), function(x){
  model<-make.bisse(tree=pruned.trees[[x]], states=epi, sampling.f=sampling.f)
  lik.4<-constrain(model, q10 ~ 0)
  start.values<-starting.point.bisse(pruned.trees[[x]])
  res.4<-find.mle(lik.4, start.values[argnames(lik.4)])
  cat("job #", x, "finished\n", sep="")
  return(res.4)
})

fit.4.est.param.table<-as.data.frame(do.call("rbind", 
      lapply(fit.4, function(x)( c(x$par, npar=length(x$par), lnl=x$lnLik)))))

#write the table to a file
output<-paste(Sys.getenv("HOME"), "/fit.4.txt", sep='')
write.table(fit.4.est.param.table, file=output, row.names=F, quote=F, sep="\t")

#Mk2 model
fit.5<-lapply(1:length(pruned.trees), function(x){
  model<-make.bisse(tree=pruned.trees[[x]], states=epi, sampling.f=sampling.f)
  lik.5<-constrain(model, lambda0 ~ lambda1, mu0 ~ mu1)
  start.values<-starting.point.bisse(pruned.trees[[x]])
  res.5<-find.mle(lik.5, start.values[argnames(lik.5)])
  cat("job #", x, "finished\n", sep="")
  return(res.5)
})

fit.5.est.param.table<-as.data.frame(do.call("rbind", 
      lapply(fit.5, function(x)( c(x$par, npar=length(x$par), lnl=x$lnLik)))))

#write the table to a file
output<-paste(Sys.getenv("HOME"), "/fit.5.txt", sep='')
write.table(fit.5.est.param.table, file=output, row.names=F, quote=F, sep="\t")

#Mk1 model
fit.6<-lapply(1:length(pruned.trees), function(x){
  model<-make.bisse(tree=pruned.trees[[x]], states=epi, sampling.f=sampling.f)
  lik.6<-constrain(model, lambda0 ~ lambda1, mu0 ~ mu1, q01 ~ q10)
  start.values<-starting.point.bisse(pruned.trees[[x]])
  res.6<-find.mle(lik.6, start.values[argnames(lik.6)])
  cat("job #", x, "finished\n", sep="")
  return(res.6)
})

fit.6.est.param.table<-as.data.frame(do.call("rbind", 
      lapply(fit.6, function(x)( c(x$par, npar=length(x$par), lnl=x$lnLik)))))

#write the table to a file
output<-paste(Sys.getenv("HOME"), "/fit.6.txt", sep='')
write.table(fit.6.est.param.table, file=output, row.names=F, quote=F, sep="\t")

#no state-dependent diversification, no reversals
fit.7<-lapply(1:length(pruned.trees), function(x){
  model<-make.bisse(tree=pruned.trees[[x]], states=epi, sampling.f=sampling.f)
  lik.7<-constrain(model, lambda0 ~ lambda1, mu0 ~ mu1, q10 ~ 0)
  start.values<-starting.point.bisse(pruned.trees[[x]])
  res.7<-find.mle(lik.7, start.values[argnames(lik.7)])
  cat("job #", x, "finished\n", sep="")
  return(res.7)
})

fit.7.est.param.table<-as.data.frame(do.call("rbind", 
     lapply(fit.7, function(x)( c(x$par, npar=length(x$par), lnl=x$lnLik)))))

#write the table to a file
output<-paste(Sys.getenv("HOME"), "/fit.7.txt", sep='')
write.table(fit.7.est.param.table, file=output, row.names=F, quote=F, sep="\t")
