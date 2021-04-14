#' Maximum Likelihood Estimate for psi
#'
#'
#' Numerically searches for the MLE of psi as the root of equation
#' K<-sum(psi/psi+i-1) for i<-1:n, where K is the observed number of
#' different species in the sample. An accepted psi sets value of the right side
#' of the equation within R's smallest possible value of the actual value of K.
#' @param abund An abundance vector.
#' @keywords maximum likelihood estimate psi
#' @details Numerically searches for the MLE of psi as the root of equation
#' K<-sum(psi/psi+i-1) for i<-1:n, where K is the observed number of
#' different species in the sample. An accepted psi sets value of the right side
#' of the equation within R's smallest possible value of the actual value of K.
#' @return psi   Maximum Likelihood Estimate of psi
#' @return Asymptotic confidence interval
#' @export
#' @examples
#' MLEp(abundances(table(c(1,2,2))))
#'
#' set.seed(1000)
#' x<-rPD(10000,2000)
#' MLEp(abundances(table(x)))



MLEp<- function(abund) {
  n<-sum(as.integer(names(abund))*abund)
  k<-sum(abund)
  psi<-1
  asum<-0
  last<-psi/2

  while(abs(asum-k)>.Machine$double.xmin) {


    if (asum<k && last==psi/2) {
      last<-psi
      psi<-2*psi
    } else if (asum<k){
      psi1<-psi
      psi<- psi + abs(psi-last)/2
      last<-psi1
    } else if (asum>k && last==psi/2){
      last<-psi
      psi<-0.5*psi
    } else {
      psi1<-psi
      psi<- psi - abs(psi-last)/2
      last<-psi1
    }

    if(last==psi) {break}

    asum<-sum(psi/(psi+seq(0,n-1)))
  }
  sum1=sum(1/(psi + seq(0,n-1)))
  sum2=sum(1/(psi + seq(0,n-1))**2)
  lower.bound<-psi + qnorm(0.025) / sqrt(n*(sum1/psi - sum2))
  upper.bound<-psi + qnorm(0.975) / sqrt(n*(sum1/psi - sum2))
  return(list("psi"=psi,"Asymptotic confidence interval"=c(lower.bound, upper.bound)))
}

