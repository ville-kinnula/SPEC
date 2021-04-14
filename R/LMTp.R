#' Lagrange Multiplier Test for psi
#'
#' Performs the Lagrange Multiplier test for an abundance vector under partition exchangeability.
#' @param psi Target psi to be tested. psi = "a" for absolute value 1, "r" for relative value n (sample size), or any positive number.
#' @param abund An abundance vector.
#' @return A p-value.
#' @keywords score test
#' @details
#' U(psi_0)^2/I(psi_0),
#' where U is the log-likelihood function of psi and I is its Fisher information.
#' The statistic follows khi-squared distribution with 1 degree of freedom
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
  return(1-pchisq((k/psi - asum)**2/(asum/psi - sum(1/(psi + seq(0,n-1))**2)),1))
}
