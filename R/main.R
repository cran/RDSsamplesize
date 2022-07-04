#' Calculate the  accumulated sample size distribution by each wave.
#'
#' @param n0 Number of seeds.
#' @param m Number of coupons issued to each participant.
#' @param maxT Planned field period.
#' @param p_list_vec A vector of recruitment rates.
#' @param tol Accuracy loss limit.
#' @return a list consisting of the following elements:
#'\item{P_tau_list}{vector; a vector of extinction probability at each wave.}
#'\item{Fk}{list; probability mass function of the accumulated sample size by each wave, k=1,...,maxT.}
#' @useDynLib RDSsamplesize, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @examples
#' result <- calSize(n0=10,m=3,maxT=9,p_list_vec=rep(0.3,9),tol=0.005)
#'
#' @export
calSize<-function(n0,m,maxT,p_list_vec,tol){
  size(n0,m,maxT,p_list_vec,tol)
}
