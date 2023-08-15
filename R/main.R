#' Calculating the accumulated sample size distribution by each wave.
#'
#' @param s scalar; Number of seeds to initiate the sampling process.
#' @param c scalar; Number of coupons issued to each participant.
#' @param maxWave scalar; Planned field period scaled by wave, which does not include the initial round of recruiting seeds.
#' @param rr scalar or vector; a (constant) recruitment rate or a vector of length \emph{maxWave}, listing varying recruitment rates at each wave.
#' The recruitment rate represents the average coupon use rate. For example, if \emph{rr} is a vector,
#' the \emph{w}th element is the ratio of the number of successful recruits brought into the study at wave \emph{w} by their recruiters (participants from wave \emph{w-1})
#' to the total number of coupons issued to those recruiters, where \emph{w} ranges from 1 to \emph{maxWave}. Seeds are counted as participants at Wave 0.
#' @param bruteMC logical; If TRUE then use a brute force Monte Carlo approach to obtain empirical data and estimate sample size distribution;
#' If FALSE then compute the theoretical results of sample size distribution using an approximation algorithm.
#' @param tol scalar; Accuracy loss limit control, which is set up for the approximation algorithm when \emph{bruteMC}=FALSE, with default of 0.025.
#' This parameter determines the acceptable level of accuracy loss in the approximate computation of the sample size distribution.
#' @return a list consisting of the following elements:
#'\item{Pr_Extinction_list}{vector; a vector of extinction probabilities, i.e., probability of not recruiting any new participants at each wave.}
#'\item{Pr_Size_by_Wave_w}{list; probability mass function and complementary cumulative distribution function of
#'attaining a certain sample size (including seeds) by each wave, w=1,...,maxWave. The round of seed collection is counted as wave 0.}
#' @useDynLib RDSsamplesize, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats rbinom ecdf
#' @references
#'
#' Raychaudhuri, Samik. \emph{Introduction to monte carlo simulation}, 2008 Winter simulation conference. IEEE, 2008.
#' @examples
#' x <- calSize(s=10,c=3,maxWave=9,rr=0.3,bruteMC=FALSE,tol=0.025)
#'
#' @export
calSize<-function(s,c,maxWave,rr,bruteMC,tol=0.025){
  if(bruteMC==FALSE){
    x<-size(s,c,maxWave,rr,tol)
    for(i in 2:length(x)){
      tb<-x[[i]]
      tb[,1]<-tb[,1]+s#add seeds
      #calculate CCDF
      ccdf<-1-cumsum(tb[,2])
      tb<-cbind(tb,c(1,ccdf[1:(length(ccdf)-1)]))
      colnames(tb)<-c("Accumulated size","PMF","CCDF")
      x[[i]]<-tb
    }
    names(x)<-c("Pr_Extinction",paste0("Pr_Size_by_Wave_",1:(length(x)-1)))
  }
  else{
    if(length(rr)==1&maxWave>1)
      rr=rep(rr,maxWave)
    sims <- 100000
    ss<-matrix(NA,nrow=sims,ncol=maxWave+1)
    ss[,1]<-s
    for(wave in 1:maxWave){
      ss[,wave+1]<-rbinom(sims,ss[,wave]*c,rr[wave])
    }
    cumss<-t(apply(ss,1,function(x)cumsum(x)))
    x<-list()
    l<-apply(ss,2,function(x)mean(x==0))[2:ncol(ss)]
    x[[1]]<-l-c(0,l[1:(length(l)-1)])
    for(wave in 1:maxWave){
      f<-ecdf(cumss[,wave+1])
      CDF=f(s:(s*c^wave))
      tb<-data.frame(`Accumulated size`=s:(s*c^wave),PMF=CDF-c(0,CDF[1:(length(CDF)-1)]),
                     CCDF=1-CDF)
      tb<-tb[tb$CCDF>0,]
      colnames(tb)[1]<-"Accumulated size"
      x[[wave+1]]<-tb
    }
    names(x)<-c("Pr_Extinction",paste0("Pr_Size_by_Wave_",1:(length(x)-1)))
  }
  return(x)
}


#' Summarizing the sample size estimation.
#'
#' @param x an object class of "RDSsamplesize", results of estimated sample size distribution of a call to 'calSize'.
#' @param n integer; target sample size.
#' @return a table presenting the probability of the accumulated sample size (including seeds) reaching at least \emph{n}
#' by each wave, w=1,..., \emph{maxWave}
#' @useDynLib RDSsamplesize, .registration = TRUE
#' @examples
#' x <- calSize(s=10,c=3,maxWave=9,rr=0.3,bruteMC=FALSE,tol=0.025)
#' nprobw(x,n=100)
#' @export
nprobw<-function(x,n){
  tb<-matrix(NA,nrow = 1, ncol=length(x)-1)
  for(i in 2:length(x)){
    ccdf<-x[[i]]
    if(!n%in%ccdf[,1])
      tb[i-1]<-0
    else
      tb[i-1]<-ccdf[match(n,ccdf[,1]),3]
  }
  if(length(x)>1)
    colnames(tb)<-c("By wave 1", paste0("wave ",2:(length(x)-1)))
  else
    colnames(tb)<-"By wave 1"
  rownames(tb)<-paste0("Pr(size>=",n,")")
  return(tb)
}

