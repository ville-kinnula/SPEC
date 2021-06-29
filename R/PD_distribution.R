#' Poisson Dirichlet distribution
#'
#' Probability of a data vector x given by Ewens' sampling formula. The higher the
#' dispersal parameter the psi is, the higher the amount of distinct observed
#' species will be. In terms of the paintbox process, a high psi increases the
#' size of the continuous part p_0 of the process, while a low psi will increase
#' the size of the discrete parts p_>0.

#' @usage dPD(abund, psi)
#' @param abund An abundance vector.
#' @param psi Dispersal parameter. Accepted values are positive numbers, "a" for absolute value psi=1 by default, or "r" for relative value psi equals sample size.
#' @return dPD returns the probability of the abundance vector of the data vector x, given dispersal parameter psi.
#' @keywords Poisson-Dirichlet distribution
#' @export
#' @references W.J. Ewens, The sampling theory of selectively neutral alleles, Theoretical Population Biology, Volume 3, Issue 1,
#' 1972, Pages 87-112, ISSN 0040-5809, https://doi.org/10.1016/0040-5809(72)90035-4.
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





#' Sampling from he Dirichlet-Poisson Distribution
#'
#' rPD samples from the PD distribution by simulating the Hoppe urn model.
#' @usage rPD(n, psi)
#' @param n number of observations.
#' @param psi dispersal parameter.
#' @return rPD returns a list with a sample of size n from the Hoppe urn model with parameter psi, along with its table of frequencies.
#' given parameter psi
#' @keywords Poisson-Dirichlet distribution
#' @export
#' @references Hoppe, F.M. The sampling theory of neutral alleles and an urn model in population genetics.
#'  J. Math. Biology 25, 123–159 (1987). https://doi.org/10.1007/BF00276386
#' @examples
#' set.seed(111)
#' s <- rPD(100,5)

rPD <- function(n, psi) {
  PDsample <- c(1)
  thesample<-c(1)
  for (i in 2:n) {
    newInt <- sample.int(length(PDsample)+1, 1, prob=c(PDsample/(i-1+psi), psi/(i-1+psi)))
    if (newInt==length(PDsample)+1) { PDsample<-append(PDsample,1) } else { PDsample[newInt]<-PDsample[newInt]+1}
    thesample<-append(thesample, newInt)
  }

  return(thesample)
}
