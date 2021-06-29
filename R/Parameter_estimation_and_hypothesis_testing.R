#' Maximum Likelihood Estimate for psi
#'
#'
#' Numerically searches for the MLE of psi as the root of equation
#' K=sum(psi/psi+i-1) for i<-1:n, where K is the observed number of
#' different species in the sample. An accepted psi sets the value of the right side
#' of the equation within R's smallest possible value of the actual value of K.
#' Returns a list that contains the estimate "psi" and "Asymptotic confidence interval".
#' The confidence interval is based on the asymptotic distribution of maximum likelihood estimators.
#' @param abund An abundance vector.
#' @keywords maximum likelihood estimate psi
#' @details Numerically searches for the MLE of psi as the root of equation
#' K<-sum(psi/psi+i-1) for i<-1:n, where K is the observed number of
#' different species in the sample. An accepted psi sets value of the right side
#' of the equation within R's smallest possible value of the actual value of K.
#' @return Returns a list containing the MLE in $psi and an asymptotic confidence interval in
#' $´Asymptotic confidence interval´
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
  lower.bound<-psi + stats::qnorm(0.025) / sqrt(n*(sum1/psi - sum2))
  upper.bound<-psi + stats::qnorm(0.975) / sqrt(n*(sum1/psi - sum2))
  return(list("psi"=psi,"Asymptotic confidence interval"=c(lower.bound, upper.bound)))
}






#' Bootstrap confidence interval for the MLE of psi
#'
#' A bootstrapped confidence interval for the Maximum Likelihood Estimate for psi.
#' @param x A data vector.
#' @export
#' @return Lower and upper bounds of the confidence interval.


MLEp.bsci<-function(x) {
  n_bootstrap<-1000
  n<-floor(0.8*length(x))
  bootsrap_psis<-c()
  for (i in 1:n_bootstrap) {
    bootsrap_psis<-append(bootsrap_psis, MLEp(table(table(sample(x,n))))$psi)
  }
  return(c(MLEp(table(table(x)))$psi, stats::quantile(bootsrap_psis,c(0.025,0.975))))
}



#' Lagrange Multiplier Test for psi
#'
#' Performs the Lagrange Multiplier test for an abundance vector under partition
#' exchangeability. Returns a p-value for the hypothesis that the input data
#' vector stems from a population with the input dispersal parameter.
#' @param psi Target psi to be tested. psi = "a" for absolute value 1, "r" for relative value n (sample size), or any positive number.
#' @param abund An abundance vector.
#' @return A p-value.
#' @keywords score test
#' @details
#' U(psi_0)^2/I(psi_0),
#' where U is the log-likelihood function of psi and I is its Fisher information.
#' The statistic follows chi-squared distribution with 1 degree of freedom
#' when the null hypothesis H_0:psi<-psi_0 is true.
#' @export
#' @examples
#' set.seed(10000)
#' x<-rPD(1000, 10)
#' abund=abundances(table(x))
#' LMTp(abund, 10)
#' LMTp(abund, 15)
#' LMTp(abund, 5)
#' LMTp(abund)       #test for psi=1
#' LMTp(abund, "r")  #test for psi=n

LMTp <- function(abund, psi="a") {
  n<-sum(as.integer(names(abund))*abund)
  k<-sum(abund)
  if (psi=="a") {
    psi<-1
  } else if (psi=="r") {
    psi<-n
  } else if (!is.numeric(psi) || psi<=0) {
    return("Psi must be a positive number.")
  }
  asum<-sum(1/(psi + seq(0,n-1)))
  return(1-stats::pchisq((k/psi - asum)**2/(asum/psi - sum(1/(psi + seq(0,n-1))**2)),1))
}
