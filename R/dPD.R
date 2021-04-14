#' Poisson Dirichlet distribution
#'
#' Probability of a data vector x from Ewens' sampling formula. The higher the
#' dispersal parameter the psi is, the higher the amount of distinct observed
#' species will be. In terms of the paintbox process, a high psi increases the
#' size of the continuous part p_0 of the process, while a low psi will increase
#' the size of the discrete parts p_>0.
#'
#' @usage dPD(abund, psi)
#' @param abund An abundance vector.
#' @param psi Dispersal parameter. Accepted values are positive numbers, "a" for absolute value psi=1 by default, or "r" for relative value psi equals sample size.
#' @return dPD returns the probability of the abundance vector of the data vector x,
#' @keywords Poisson-Dirichlet distribution
#' @export
#' @examples
#' set.seed(111)
#' s <- rPD(100,5)
#' a=table(table(s))
#' dPD(a, 5)
#'
dPD <- function(abund, psi="a") {
  n<-sum(as.integer(names(abund))*abund)
  if (psi=="a") {
    psi<-1
  } else if (psi=="r") {
    psi<-n
  } else if (!is.numeric(psi) || psi<=0) {
    return("Psi must be a positive number.")
  }
  return(exp(dlPD(abund, psi)))
}
