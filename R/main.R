#' Calculate the  accumulated sample size distribution by each wave.
#'
#' @param n0 Number of seeds.
#' @param m Number of coupons issued to each participant.
#' @param maxT Planned field period.
#' @param p_list_vec A vector of recruitment rates by each wave.
#' @param tol Accuracy loss limit.
#' @return a list consisting of the following elements:
#'\item{Pr_Extinction_list}{vector; a vector of extinction probability, i.e., probability of recruiting none at each wave.}
#'\item{Pr_Size_by_Wave_t}{list; probability mass function of the accumulated sample size (excluding seeds) by each wave, t=1,...,maxT.}
#' @useDynLib RDSsamplesize, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @examples
#' result <- calSize(n0=10,m=3,maxT=9,p_list_vec=rep(0.3,9),tol=0.005)
#'
#' @export
calSize<-function(n0,m,maxT,p_list_vec,tol){
  output<-size(n0,m,maxT,p_list_vec,tol)
  names(output)<-c("Pr_Extinction",paste0("Pr_Size_by_Wave_",1:(length(output)-1)))
  return(output)
}
