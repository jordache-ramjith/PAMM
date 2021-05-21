library(ggplot2)
library(pammtools)
#remotes::install_github("csdaw/ggprism")
library(ggprism)
library(data.table)
library(dplyr)
library(survival)
library(mgcv)
library(coxme)


#################

simulation=function(ns=10,ne=200, sample=100, var=0,beta=-1,maxtime=25){sim=list()
rows=ne*sample
library(dplyr)
library(survival)
for (j in 1:ns){
  n=c(rep(seq(1:sample),each=ne))
  sim[[j]]=data.frame(n)
  
  
  set.seed(j)
  sim[[j]] = sim[[j]] %>% mutate(event=c(rep(seq(1:ne),sample)),var=var,x=c(rep(rbinom(sample,1,0.5),each=ne)),v=c(rep(rnorm(sample,0,sqrt(var)),each=ne)),beta=c(rep(beta,rows)),u2=c(rep(runif(ne*sample,0,1))))
  
  sim[[j]]$start=0
  
  for (i in 1:rows){
    sim[[j]]$time[i]=if (sim[[j]]$event[i]==1) (((-log(sim[[j]]$u2[i])/((2*exp(sim[[j]]$beta[i]*sim[[j]]$x[i]+sim[[j]]$v[i]))))))^(1/0.5) else (((-log(sim[[j]]$u2[i])/((2*exp(sim[[j]]$beta[i]*sim[[j]]$x[i]+sim[[j]]$v[i])))))+(sim[[j]]$time[i-1]^(0.5)))^(1/0.5)
  }
  
  sim[[j]] = sim[[j]] %>% group_by(n) %>% mutate(start = data.table::shift(time)) %>% ungroup()  %>% mutate(start=ifelse(is.na(start)==T,0,start),maxtime=maxtime,status=ifelse(time<maxtime,1,0),time=ifelse(status==1,time,maxtime),keep=ifelse(time-start<=0.01,"No",""))
  
  sim[[j]]=as.data.frame(sim[[j]][!sim[[j]]$keep == "No", ] %>% group_by(n) %>% mutate(event=row_number()) %>% ungroup())
  
}


return(sim)

}
#######ANALYSIS##########


simdat=list()

#no heterogeneity
simdat[[1]]=simulation(ns=200,sample=200)
#low heterogeneity
simdat[[2]]=simulation(ns=200,var=0.2,sample=200)
#moderate heterogeneity
simdat[[3]]=simulation(ns=200,var=0.5,sample=200)
#high heterogeneity
simdat[[4]]=simulation(ns=200,var=1,sample=200)
#no heterogeneity
simdat[[5]]=simulation(ns=200,maxtime = 100,sample=200)
#low heterogeneity
simdat[[6]]=simulation(ns=200,var=0.2,maxtime = 100,sample=200)
#moderate heterogeneity
simdat[[7]]=simulation(ns=200,var=0.5,maxtime = 100,sample=200)
#high heterogeneity
simdat[[8]]=simulation(ns=200,var=1,maxtime = 100,sample=200)

sumey=lapply(simdat,FUN=function(x){lapply(x, FUN=function(y){sum(y$status[y$event==200])})})
lapply(sumey, FUN=function(x){Reduce("+",x)})


dat=list()

for (j in 1:length(simdat)){
  dat[[j]]=data.frame(coef1=NA,se1=NA,sigma1=NA,pamcoef=NA,pamse=NA,pamsigma=NA,pamsigmalci=NA,pamsigmauci=NA,var=NA,coverage_cox_coef=NA,coverage_pam_coef=NA)
  for(i in 1:length(simdat[[j]])){
    dat[[j]][i,]$var=unique(simdat[[j]][[i]]$var)
    coxmod1 = coxme(Surv(start,time,status)~x+(1|n),data=simdat[[j]][[i]])
    dat[[j]][i,]$coef1=coef(coxmod1)
    dat[[j]][i,]$se1=sqrt(vcov(coxmod1))
    dat[[j]][i,]$sigma1=coxmod1$vcoef$n
    
    ped <- as_ped(
      formula     = Surv(start,time,status) ~ x,
      id          = "n",
      data        = simdat[[j]][[i]],
      transition  = "event",
      timescale   = "calendar",
      max_time = max(simdat[[j]][[i]]$time)) %>%
      mutate(
        event = as.factor(event),
        n = as.factor(n))
    
    pam<- bam(ped_status ~ x+s(tend)+s(n,bs="re"),data= ped, family = poisson(), offset = offset, discrete = TRUE)
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
    print(paste("j=",j,"i=" , i))
  }
}


sumstats=lapply(simdat[[1]],FUN=function(x){psych::describeBy(x)})


betadens=lapply(dat,FUN=function(x){ggplot(data=x)+
    geom_density(aes(x=coef1,fill="CPH frailty model"),alpha=0.5)+
    xlab("")+ylab("")+
    scale_x_continuous(breaks = seq(-2,0,0.25),limits=c(-2,0),expand=c(0,0))+
    scale_y_continuous(breaks = seq(0,8,1),limits=c(0,8.5),expand=c(0,0))+
    geom_density(aes(x=pamcoef,fill="PAMM"),alpha=0.5)+
    guides(fill=guide_legend(title=element_blank(),nrow = 1))+
    theme(legend.position = "top")+
    scale_fill_manual(values=c("#CC0000","#003366"))+
    theme_prism()+
    geom_vline(xintercept = -1,linetype="dashed",size=1.2)})

library(ggpubr)
betaplot=ggarrange(
  plotlist=betadens, labels = c("Scenario 1", "Scenario 2","Scenario 3","Scenario 4","Scenario 5", "Scenario 6","Scenario 7","Scenario 8"),nrow=4,ncol=2,vjust=0,hjust=-2.3,
  common.legend = TRUE
)

betaplot=annotate_figure(betaplot,left=text_grob("Density",face="bold",size=14,rot=90),bottom=text_grob(bquote(bold("Estimated coefficient "*hat(beta)*"")),face="bold",size=14))

betaplot

ggsave("betaplot_calendar_largeN.pdf",betaplot,width = 10, height = 7, dpi = 300,device='pdf')


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


library(ggpubr)
sigmaplot=ggarrange(
  plotlist=sigmadens, labels = c("Scenario 1", "Scenario 2","Scenario 3","Scenario 4","Scenario 5", "Scenario 6","Scenario 7","Scenario 8"),nrow=4,ncol=2,vjust=0,hjust=-2.3,
  common.legend = TRUE
)

sigmaplot=annotate_figure(sigmaplot,left=text_grob("Density",face="bold",size=14,rot=90),bottom=text_grob(bquote(bold("Estimated frailty variance "*hat(sigma)^{2}*"")),face="bold",size=14))

sigmaplot

ggsave("sigmaplot_calendar_largeN.pdf",sigmaplot,width = 10, height = 7, dpi = 300,device='pdf')

shortdat=lapply(dat, FUN=function(x){x%>% mutate(diff=coef1-pamcoef) %>%select(coef1,pamcoef,se1,pamse,sigma1,pamsigma,coverage_cox_coef,coverage_pam_coef,diff)})

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


save.image("WS_box_steffensmeier_de_boef_ct_largeN.RData")

