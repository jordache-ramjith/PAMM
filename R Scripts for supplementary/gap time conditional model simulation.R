library(ggplot2)
library(pammtools)
library(mgcv)
#remotes::install_github("csdaw/ggprism")
library(ggprism)
library(coxme)

#j maxed at 4 for baseline to stop it from getting too small

simulation=function(ns=10,ne=40, sample=100, var=0,beta=-1, lambda=1,maxtime=5,dep=FALSE){sim=list()
rows=ne*sample
library(dplyr)
library(survival)
for (j in 1:ns){
  n=c(rep(seq(1:sample),each=ne))
  sim[[j]]=data.frame(n)
  
  
  
  set.seed(j)
  sim[[j]] = sim[[j]] %>% mutate(event=c(rep(seq(1:ne),sample)),var=var,x=c(rep(rbinom(sample,1,0.5),each=ne)),v=c(rep(rnorm(sample,0,sqrt(var)),each=ne)),beta=c(rep(beta,rows)),u2=c(rep(runif(ne*sample,0,1))),lambda=lambda,dep=dep,lambda=ifelse(dep==TRUE,apply(cbind(lambda*event,lambda*4), 1, FUN=min),lambda),Ti=(-3*log(u2)/((lambda*exp(beta*x+v))))^(1/3),time=ave(Ti, n, FUN=cumsum),start=time-Ti,maxtime=maxtime,status=ifelse(time<maxtime,1,0),time=ifelse(status==1,time,maxtime),keep=ifelse(start>=time,"No",""),gap=time-start)
  

  sim[[j]]=sim[[j]][!sim[[j]]$keep == "No", ]
  
}


return(sim)

}


#######ANALYSIS##########


simdat=list()

#no dependence + no heterogeneity
simdat[[1]]=simulation(ns=500)
#no dependence plus low heterogeneity
simdat[[2]]=simulation(ns=500,var=0.2)
#no dependence plus moderate heterogeneity
simdat[[3]]=simulation(ns=500,var=0.5)
#no dependence + high heterogeneity
simdat[[4]]=simulation(ns=500,var=1)

#dependence + no heterogeneity
simdat[[5]]=simulation(ns=500,dep=TRUE,ne=50,maxtime=10)
#dependence plus low heterogeneity
simdat[[6]]=simulation(ns=500,var=0.2,dep=TRUE,ne=50,maxtime=10)
#dependence plus moderate heterogeneity
simdat[[7]]=simulation(ns=500,var=0.5,dep=TRUE,ne=50,maxtime=10)
#dependence + high heterogeneity
simdat[[8]]=simulation(ns=500,var=1,dep=TRUE,ne=50,maxtime=10)

sumey=lapply(simdat,FUN=function(x){lapply(x, FUN=function(y){sum(y$status[y$event==50])})})
lapply(sumey, FUN=function(x){Reduce("+",x)})

dat=list()

for (j in 1:length(simdat)){
  dat[[j]]=data.frame(coef1=NA,se1=NA,sigma1=NA,pamcoef=NA,pamse=NA,pamsigma=NA,pamsigmalci=NA,pamsigmauci=NA,var=NA,coverage_cox_coef=NA,coverage_pam_coef=NA)
for(i in 1:length(simdat[[j]])){
  dat[[j]][i,]$var=unique(simdat[[j]][[i]]$var)
  coxmod1 = coxme(Surv(gap,status)~x+strata(event)+(1|n),data=simdat[[j]][[i]])
  dat[[j]][i,]$coef1=coef(coxmod1)
  dat[[j]][i,]$se1=sqrt(vcov(coxmod1))
  dat[[j]][i,]$sigma1=coxmod1$vcoef$n
  
  ped <- as_ped(
    formula     = Surv(start,time,status) ~ x,
    id          = "n",
    data        = simdat[[j]][[i]],
    transition  = "event",
    timescale   = "gap",
    max_time = max(simdat[[j]][[i]]$gap)) %>%
    mutate(
      event = as.factor(event),
      n = as.factor(n))
  
  pam<- bam(ped_status ~ x+event+s(tend,by=factor(event))+s(n,bs="re"),data= ped, family = poisson(), offset = offset, discrete = TRUE)
  sm=summary(pam)
  t=gam.vcomp(pam)

  dat[[j]][i,]$pamcoef=pam$coefficients[2]
  dat[[j]][i,]$pamse=sm$se[2]
  dat[[j]][i,]$pamsigma=(t[nrow(t),1])^2
  dat[[j]][i,]$pamsigmalci=(t[nrow(t),2])^2
  dat[[j]][i,]$pamsigmauci=(t[nrow(t),3])^2
    dat[[j]][i,]$coverage_cox_coef=if(-1<=dat[[j]][i,]$coef1+1.96*dat[[j]][i,]$se1 & -1>=dat[[j]][i,]$coef1-1.96*dat[[j]][i,]$se1) 1 else 0
  dat[[j]][i,]$coverage_pam_coef=if(-1<=dat[[j]][i,]$pamcoef+1.96*dat[[j]][i,]$pamse & -1>=dat[[j]][i,]$pamcoef-1.96*dat[[j]][i,]$pamse) 1 else 0
  rm(pam,ped,sm,t)
  print(paste("j=",j,"i=",i))
}
}


sumstats=lapply(dat,FUN=function(x){psych::describeBy(x)})

betadens=lapply(dat,FUN=function(x){ggplot(data=x)+
    geom_density(aes(x=coef1,fill="CPH frailty model"),alpha=0.5)+
    xlab("")+ylab("")+
    scale_x_continuous(breaks = seq(-2,0,0.25),limits=c(-2,0),expand=c(0,0))+
    scale_y_continuous(breaks = seq(0,5,1),limits=c(0,5),expand=c(0,0))+
    geom_density(aes(x=pamcoef,fill="PAMM"),alpha=0.5)+
    guides(fill=guide_legend(title=element_blank(),nrow = 1))+
    theme(legend.position = "top")+
    scale_fill_manual(values=c("#CC0000","#003366"))+
    theme_prism()+
    geom_vline(xintercept = -1,linetype="dashed",size=1.2)})

library(ggpubr)
betaplot=ggarrange(
  plotlist=betadens, labels = c("Scenario 1", "Scenario 2","Scenario 3","Scenario 4","Scenario 5","Scenario 6","Scenario 7","Scenario 8"),nrow=4,ncol=2,vjust=0,hjust=-2.3,
  common.legend = TRUE
)

betaplot=annotate_figure(betaplot,left=text_grob("Density",face="bold",size=14,rot=90),bottom=text_grob(bquote(bold("Estimated coefficient "*hat(beta)*"")),face="bold",size=14))

ggsave("betaplot.pdf",betaplot,width = 10, height = 7, dpi = 300,device='pdf')


sigmadens=lapply(dat,FUN=function(x){ggplot(data=x)+
    geom_density(aes(x=sigma1,fill="CPH frailty model"),alpha=0.5)+
    xlab("")+ylab("")+
    scale_x_continuous(breaks = seq(0,2,0.2),limits=c(0,2),expand=c(0,0))+
    scale_y_continuous(,expand=c(0,0))+
    geom_density(aes(x=pamsigma,fill="PAMM"),alpha=0.5)+
    guides(fill=guide_legend(title=element_blank(),nrow = 1))+
    theme(legend.position = "top")+
    scale_fill_manual(values=c("#CC0000","#003366"))+
    theme_prism()+
    geom_vline(xintercept = x$var,linetype="dashed",size=1.2)})

sigmadens[[4]]

library(ggpubr)
sigmaplot=ggarrange(
  plotlist=sigmadens, labels = c("Scenario 1", "Scenario 2","Scenario 3","Scenario 4","Scenario 5","Scenario 6","Scenario 7","Scenario 8"),nrow=4,ncol=2,vjust=0,hjust=-2.3,
  common.legend = TRUE
)

sigmaplot=annotate_figure(sigmaplot,left=text_grob("Density",face="bold",size=14,rot=90),bottom=text_grob(bquote(bold("Estimated frailty variance "*hat(sigma)^{2}*"")),face="bold",size=14))

sigmaplot

ggsave("sigmaplot.pdf",sigmaplot,width = 10, height = 7, dpi = 300,device='pdf')

shortdat=lapply(dat, FUN=function(x){x%>% mutate(diff=coef1-pamcoef)%>%select(coef1,pamcoef,se1,pamse,sigma1,pamsigma,coverage_cox_coef,coverage_pam_coef,diff)})

library(data.table)
longshortdat<-rbindlist(lapply(1:length(shortdat)
                     , function(x){ setDT(shortdat[[x]])[
                       , id:=x]})
              , use.names=TRUE, fill=TRUE)

longshortdat$id=as.factor(longshortdat$id)

table1=longshortdat %>% group_by(id) %>% summarise(coefcox=mean(coef1),SEM1=sd(coef1),cov1=mean(coverage_cox_coef)*100,sigma1=median(sigma1),coef2=mean(pamcoef),SEM2=sd(pamcoef),cov2=mean(coverage_pam_coef)*100,sigma2=median(pamsigma),mdiff=mean(diff),sediff=sd(diff))

table1

library(xtable)
print(xtable(table1,digits=c(0,0,2,2,1,2,2,2,1,2,1,1), display=c("f","f", "f","f","f","f", "f","f","f","f","E","E")),floating=FALSE,latex.environments=NULL,booktabs=TRUE, include.rownames=FALSE, math.style.exponents = TRUE)


save.image("WS_box_steffensmeier_de_boef.RData")




