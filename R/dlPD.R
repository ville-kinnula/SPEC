#' Poisson Dirichlet distribution
#'
#' LogProbability of a data vector x from Ewens' sampling formula is given below.
#' accepts either a raw data vector x or its frequency vector, table(x). The higher the
#' dispersal parameter the psi is, the higher the amount of distinct observed
#' species will be. In terms of the paintbox process, a high psi increases the
#' size of the continuous part p_0 of the process, while a low psi will increase
#' the size of the discrete parts p_1, ... p_k.
#'
#' @usage dlPD(abund, psi)
#' @param abund An abundance vector.
#' @param psi Dispersal parameter. Accepted values are positive numbers, "a" for absolute value psi=1 by default, or "r" for relative value psi equals oriignal sample size.
#' @return dlPD returns the log-probability of the abundance vector of the data vector x,
#' @keywords Poisson-Dirichlet distribution
#' @export
#' @examples
#' set.seed(111)
#' s <- rPD(100,5)
#' a=abundances(table(s))
#' dlPD(a, 5)
#'
dlPD <- function(abund, psi="a") {
  rho<-abund
  t<-as.integer(names(abund))
  n<-sum(t*rho)

  if (psi=="a") {
    psi<-1
  } else if (psi=="r") {
    psi<-n
  } else if (!is.numeric(psi) || psi<=0) {
    return("Psi must be a positive number.")
  }

  product<-sum(rho*log(psi/t) - lfactorial(rho))


  return(lfactorial(n) - lognRF(psi,n) + product)
}
