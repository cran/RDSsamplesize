## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning=FALSE,
  message=FALSE
)

## ----setup--------------------------------------------------------------------
library(RDSsamplesize)
library(microbenchmark)
library(ggplot2)
library(dplyr)
library(latex2exp)

## -----------------------------------------------------------------------------
A1 <-calSize(s=10,c=3,maxWave=9,rr=0.4,bruteMC=FALSE,tol=0.025)
A2 <-calSize(s=10,c=3,maxWave=9,rr=0.4,bruteMC=TRUE)

## -----------------------------------------------------------------------------
microbenchmark(A1 =calSize(s=10,c=3,maxWave=9,rr=0.4,bruteMC=FALSE,tol=0.025),
               A2 =calSize(s=10,c=3,maxWave=9,rr=0.4,bruteMC=TRUE),
               times=10)

## -----------------------------------------------------------------------------
B <- calSize(s=10,c=3,maxWave=3,c(0.3,0.2,0.1),bruteMC=FALSE,tol=0.025)

## -----------------------------------------------------------------------------
round(nprobw(A1,n=210),3) # recruit 200 more individuals besides 10 seeds

## -----------------------------------------------------------------------------
round(nprobw(A2,n=210),3) # recruit 200 more individuals besides 10 seeds

## -----------------------------------------------------------------------------
round(nprobw(B,n=15),3) # recruit 5 more individuals besides 10 seeds
round(nprobw(B,n=50),3) # recruit 40 more individuals besides 10 seeds

## ----fig.width=5, fig.height=3, fig.align='center'----------------------------
plot_survival_prob<-function(P_tau_list){
  P_tau<-data.frame(Wave=1:length(P_tau_list),PMF=P_tau_list)
  P_tau<-P_tau%>%mutate(CDF=cumsum(PMF),CCDF=1-CDF)
  print(ggplot(P_tau,aes(x=Wave,y=CCDF))+
    geom_point()+
    geom_line(linetype=2)+
    scale_y_continuous(breaks=seq(0,1,by=.1),limits = c(min(0.9,min(P_tau$CCDF)),1))+
    xlab('Wave, w')+
    ylab('')+
    scale_x_continuous(breaks =1:length(P_tau_list),minor_breaks = NULL)+
    theme_minimal()+
    theme(panel.grid.minor.y = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))+
    ggtitle(TeX('$Pr(\\tau > w)$')))
}

### design A (with theoretical results) 
plot_survival_prob(P_tau_list=A1[[1]])
### design A (with empirical results) 
plot_survival_prob(P_tau_list=A2[[1]])
### design B
plot_survival_prob(P_tau_list=B[[1]])


## ----fig.width=7, fig.height=5, fig.align='center'----------------------------
plot_sample_size<-function(x){
  info<-data.frame(Wave=as.character(),Zk=as.integer(),PMF=as.double(),CCDF=as.double())
  for (i in 2:length(x)){
    info<-rbind(info,cbind(Wave=i-1,x[[i]]))
  }
  colnames(info)<-c('Wave','Size','PMF','CCDF')
  info$Wave<-paste("Wave",info$Wave)
  print(ggplot(info,aes(x=Size,y=CCDF))+
    facet_wrap(~Wave,scales = "free")+
    geom_point(size=0.5)+
    xlab('n')+
    ylab('')+
    theme_minimal()+
    theme(panel.grid.minor.y = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))+
    ggtitle(TeX('Pr( Accmulated sample size (including seeds) $\\geq n | k$-th wave)')))
}

### Design A (with theoretical results) 
plot_sample_size(A1)
### Design A (with empirical results)
plot_sample_size(A2)
### Design B 
plot_sample_size(B)


