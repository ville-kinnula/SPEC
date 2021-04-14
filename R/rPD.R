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


