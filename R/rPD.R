#' Sampling from he Dirichlet-Poisson Distribution
#'
#' rPD samples from the PD distribution by simulating the Hoppe urn model. Time
#' complexity is something like O(n^2log(psi)). A sample of 10^5 samples
#' takes 15-20 seconds with a psi<100.

#' @usage rPD(n, psi)
#' @param n number of observations.
#' @param psi dispersal parameter.
#' @return rPD returns a list with a sample of size n from the Hoppe urn model with parameter psi, along with its table of frequencies.
#' given parameter psi
#' @keywords Poisson-Dirichlet distribution
#' @export
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

rPD.slow <- function(n, psi) {
  n10<-10*n
  PDsample <- c(1)
  thesample<-c(1)
  for (i in 2:n10) {
    newInt <- sample.int(length(PDsample)+1, 1, prob=c(PDsample/(i-1+psi), psi/(i-1+psi)))
    if (newInt==length(PDsample)+1) { PDsample<-append(PDsample,1) } else { PDsample[newInt]<-PDsample[newInt]+1}
    thesample<-append(thesample, newInt)
  }

  return(sample(thesample, n))
}

#' Simulated approximation of the DP(psi)
#'
#' Algorithm A from Zarepour, M., and Al Labadi, L. (2012). On a Rapid Simulation of the Dirichlet
#' Process. Statistics and Probability Letters, 82, 916-924.
#' https://www.sciencedirect.com/science/article/pii/S0167715212000302?via%3Dihub
#' @param n Number of allowed unique species.
#' @param psi Dispersal parameter.
#' @export
#' @examples
#' m<- 10000
#' psi<- 100
#' p<- cylSet(m, psi)
#' # Sums up to 1.
#' sum(p)
#'
#' # drawing a sample from PD(psi) of size 1000
#' s<- sample.int(m, size=1000, replace = TRUE, prob = p)
#'
#' # MLE based on the sample should be close to 100, the used psi.
#' MLEp(table(table(s)))
#'

cylSet<- function(n, psi) {
  y<- runif(n)
  e<- rexp(n+1)

  gam<- sapply(1:(n+1), function(x) sum(e[1:x]))
  #gamsum<-0
  #gam<-list()
  #for(ei in e) {
  #  gamsum<- gamsum + ei
  #  gam<- list(gam, list(gamsum))
  #}
  #gam<-unlist(gam)
  Gn<- sapply(1:n, function(x) qgamma(1-gam[x]/gam[n+1], psi/n))
  probs<- Gn / sum(Gn)

  return(probs)
}


#' Sampling from he Dirichlet-Poisson Distribution
#'
#' rPDa samples from the PD distribution by simulating an approximation
#' of the cylinder set probabilities of PD(psi)
#'
#' @usage rPDa(n, psi)
#' @param n number of observations.
#' @param psi dispersal parameter.
#' @return rPDa returns a list with a sample of size n from the Hoppe urn model with parameter psi, along with its table of frequencies.
#' given parameter psi
#' @keywords Poisson-Dirichlet distribution
#' @export
#' @examples
#' set.seed(111)
#' s <- rPDa(100,5)
#'

rPDa<- function(n, psi) {
  p<- cylSet(max(n*10, psi*10), psi)
  return(sample.int(length(p), replace=TRUE, prob = p))
}

