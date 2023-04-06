## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning=FALSE,
  message=FALSE
)

## ----setup--------------------------------------------------------------------
library(RDSsamplesize)
library(ggplot2)
library(dplyr)
library(latex2exp)

## -----------------------------------------------------------------------------
########## Setup ##########
n0<-10 # number of seeds
m<-3 # number of coupons issued to each recruiter
maxT<-9 # planned field period
p_list_vec<-rep(0.3,maxT) # recruitment rate at each wave
tol<-0.005 # tolerance value for accuracy loss

########## Realization ##########
results<-calSize(n0,m,maxT,p_list_vec,tol)

## ----fig.width=5, fig.height=3, fig.align='center'----------------------------
### Extinction probability
P_tau_list<-results[[1]]
P_tau<-data.frame(Wave=1:length(P_tau_list),PMF=P_tau_list)
P_tau<-P_tau%>%mutate(CDF=cumsum(PMF),CCDF=1-CDF)
ggplot(P_tau,aes(x=Wave,y=CCDF))+
  geom_point()+
  geom_line(linetype=2)+
  scale_y_continuous( breaks=seq(0,1,by=.1),limits = c(min(P_tau$CCDF),1))+
  xlab('Wave, t')+
  ylab('')+
  scale_x_continuous(breaks =1:length(P_tau_list),minor_breaks = NULL)+
  theme_minimal()+
  theme(panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle(TeX('$Pr(\\tau > t)$'))



## ----fig.width=7, fig.height=5, fig.align='center'----------------------------
info<-data.frame(Wave=as.integer(),Zk=as.integer(),PMF=as.double())
for (i in 1:maxT){
  info<-rbind(info,cbind(Wave=i,results[[i+1]]))
}
colnames(info)<-c('Wave','Size','PMF')
info<-info%>%group_by(Wave)%>%mutate(CDF=cumsum(PMF),CCDF=1-CDF)

### CCDF of accumulated sample size 
ggplot(info,aes(x=Size,y=CCDF))+
  facet_wrap(~Wave,scales = "free")+
  geom_point(size=0.5)+
  xlab('n')+
  ylab('')+
  theme_minimal()+
  theme(panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('Pr( Accmulated sample size (excluding seeds) > n | k-th wave)')





