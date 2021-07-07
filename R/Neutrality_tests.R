#' Watterson's homozygosity test
#'
#' This function performs a statistical test on the null hypothesis thata given
#' sample's underlying distribution is the Poisson-Dirichlet distribution. It
#' calculates a test statistic that is then used to gain a p-value from an
#' empirical distribution of the statistic from simulated samples from a
#' PD distribution.
#' @param x A data vector.
#' @param rounds How many samples are simulated to obtai nthe empirical distribution.
#' @returns A p-value.
#' @export
#' @references Watterson, G.A., (1978), The homozygosity test of neutrality. Genetics. 88(2):405-417.
#' @examples
#' ##Test whether a typical sample follows PD:
#' x<-rPD(1000,10)
#' WattersonsW(x, 1000)
#'
#' ##Test whether a very atypical sample where frequencies of different values
#' ## are similar:
#'
#' x<-c(rep(1, 200), rep(2, 200), rep(3, 200), rep(4, 200), rep(5,200))
#' WattersonsW(x,1000)

WattersonsW<-function(x, rounds) {
  freq<-table(x)
  n<-length(x)
  W<-sum(freq**2)/n**2
  psi<-MLEp(table(table(x)))$psi
  W.emp<- sapply(1:rounds, function(x) sum(table(rPD(n, psi))**2)/n**2 )

  return(sum(W>W.emp)/rounds)
}




FusF<- function(x, rounds) {
  n<-length(x)
  k<-length(unique(x))
  psi<-(MLEp(table(table(x)))$psi)
  S<-sum(unlist(sapply((k+1):n, function(y) abs(gmp::Stirling1(n, y)*psi**y)))/gmp::as.bigz(exp(lognRF(psi, n))))
  Fs<- Rmpfr::log1pexp(Rmpfr::.bigq2mpfr(S/(1-S)))

  Fs.emp<-c()
  for (i in 1:rounds) {
    k<-length(unique(rPD(n, psi)))
    Si<-sum(unlist(sapply((k+1):n, function(x) abs(gmp::Stirling1(n, x)*psi**x)))/gmp::as.bigz(exp(lognRF(psi, n))))
    Fs.emp<-append(Fs.emp, Rmpfr::log1pexp(Rmpfr::.bigq2mpfr(Si/(1-Si))))
  }

  return(sum(Fs>Fs.emp)/rounds)
}
