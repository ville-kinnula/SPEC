#' Maximum Likelihood Estimate for psi
#'
#'
#' Numerically searches for the MLE of psi given the frequencies of the
#' frequencies of a data vector, called an abundance vector.
#' @param abund An abundance vector.
#' @keywords maximum likelihood estimate psi
#' @details Numerically searches for the MLE of psi as the root of equation
#' \deqn{K=\sum_{i=1}^n\psi/(\psi+i-1),} where \eqn{K} is the observed number of
#' different species in the sample. The right side of the equation is strictly
#' increasing when \eqn{\psi>0}, so a binary search is used to find the root.
#' An accepted \eqn{\psi} sets value of the right side
#' of the equation within R's smallest possible value of the actual value of \eqn{K}.
#' @return Returns a list that contains the estimate "psi" and "Asymptotic confidence interval".
#' The confidence interval is based on the asymptotic distribution of maximum likelihood estimators.
#' @export
#' @references W.J. Ewens, The sampling theory of selectively neutral alleles, Theoretical Population Biology, Volume 3, Issue 1,
#' 1972, Pages 87-112, ISSN 0040-5809, <https://doi.org/10.1016/0040-5809(72)90035-4>.
#' @examples
#' ##Find the MLE of psi of the vector (1,2,2).
#' ##The frequencies of the frequencies of the data vector are given as input:
#' MLEp(table(table(c(1,2,2))))
#'
#' ##Find the MLE of psi of a sample from the Poisson-Dirichlet distribution:
#' set.seed(1000)
#' x<-rPD(n=10000, psi=100)
#' MLEp(table(table(x)))



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
#' A bootstrapped confidence interval for the Maximum Likelihood Estimate for
#' psi based on a 100 bootstrap rounds with 80% of data used at a time.
#' @param x A data vector.
#' @param level Level of confidence interval as number between 0 and 1.
#' @export
#' @return The MLE of psi as well as lower and upper bounds of the bootstrap
#' confidence interval.
#' @examples
#' ## Find a 95% -confidence interval for the MLE of psi given a sample from the
#' ## Poisson-Dirichlet distribution:
#' x<-rPD(n=10000, psi=100)
#' MLEp.bsci(x, 0.95)
#'


MLEp.bsci<-function(x, level) {
  n_bootstrap<-1000
  n<-floor(0.8*length(x))
  bootsrap_psis<-c()
  for (i in 1:n_bootstrap) {
    bootsrap_psis<-append(bootsrap_psis, MLEp(table(table(sample(x,n))))$psi)
  }
  bounds<-c((1-level)/2, 1-(1-level)/2)
  return(c(MLEp(table(table(x)))$psi, stats::quantile(bootsrap_psis,bounds)))
}



#' Lagrange Multiplier Test for psi
#'
#' Performs the Lagrange Multiplier test for an abundance vector under partition
#' exchangeability. Returns a p-value for the hypothesis that the input data
#' vector stems from a population with the input dispersal parameter.
#' @param psi Target \eqn{\psi} to be tested. `psi` = "a" for absolute value 1,
#' "r" for relative value \eqn{n} (sample size), or any positive number.
#' @param abund An abundance vector.
#' @return A p-value.
#' @references Radhakrishna Rao, C, (1948), Large sample tests of statistical
#' hypotheses concerning several parameters with applications to problems of
#' estimation. Mathematical Proceedings of the Cambridge Philosophical Society,
#'  44(1), 50-57. <https://doi.org/10.1017/S0305004100023987>
#' @keywords score test
#' @details \deqn{U(\psi_0)^2 / I(\psi_0)},
#' where \eqn{U} is the log-likelihood function of \eqn{\psi} and \eqn{I} is its Fisher information.
#' The statistic follows \eqn{\chi^2}-distribution with 1 degree of freedom
#' when the null hypothesis \eqn{H_0:\psi=\psi_0} is true.
#' @export
#' @examples
#' ## Test the psi of a sample fro mthe Poisson-Dirichlet distribution:
#' set.seed(10000)
#' x<-rPD(1000, 10)
#' ## Find the abundances of the data vector:
#' abund=abundances(table(x))
#' ## Test for the psi that was used, as well as a high one and a low one:
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
