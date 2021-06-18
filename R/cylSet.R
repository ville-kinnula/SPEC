#' Simulated approximation of rPD
#'
#' Algorithm A from Zarepour, M., and Al Labadi, L. (2012). On a Rapid Simulation of the Dirichlet
#' Process. Statistics and Probability Letters, 82, 916-924.
#' https://www.sciencedirect.com/science/article/pii/S0167715212000302?via%3Dihub
#' @param n Number of allowed unique species.
#' @param psi Dispersal parameter.
#' @export
#'
#'

cylSet<- function(n, psi) {
  y<- runif(n)
  e<- rexp(n+1)
  gam<- sapply(1:(n+1), function(x) sum(e[1:x]))
  Gn<- sapply(1:n, function(x) qgamma(1-gam[x]/gam[n+1], psi/n))
  probs<- Gn / sum(Gn)

  return(probs)
}
