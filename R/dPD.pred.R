#' Updated predictive density from Ewens sampling formula
#'
#' Calculates predictive probability conditional on training data probability.
#' x is test data, other parameters come  from training data object.
#' coeffs avoids calculating n!/Stirlings(psi, n) again for training data.
#' @export

dPD.pred <- function(x, psi, freqs0, coeffs, tr.prob) {
  abds<-abundances(table(x), freqs0)
  rho<-abds
  t<-as.integer(names(abds))
  n<-sum(rho*t)

  product<-sum(rho*log(psi/t) - lfactorial(rho))
  #coeff<- sum(log(seq(n-length(x)+1, n))) - sum(log(seq(psi + n - length(x), psi + n-1)))

  #return(coeffs + coeff + product - tr.prob)
  #product<-sum(rho*log(psi/t) - lfactorial(rho))


  return(lfactorial(n) - lognRF(psi,n) + product - tr.prob)
}
