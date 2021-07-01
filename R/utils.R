#' Vector of frequencies of frequencies
#'
#' A function to calculate frequencies of frequencies or the combined
#' frequencies of frequencies of two vectors.
#' @param freqs A frequency table of data vector `x`.
#' @param freqs0 A second optional frequency table to be merged to the `freqs` parameter.
#' @keywords abundances
#' @details  `Freqs0` is optional. It is
#' used in the case of calculating abundances of test data when the frequencies of
#' the training data are already known. `Freqs0` is `table(x0)` ,where `x0` is training
#' data. Abundances of any kind of data vector x can be calculated with
#' `table(table(x))`.
#' @return This function returns a named vector that is used to calculate
#' probabilities and make classifications.
#' @export
#' @examples
#' set.seed(111)
#' x<-rpois(10,10)
#' abundances(table(x))
#'
#' y<-rpois(2,10)
#' abundances(table(x), table(y))

abundances<-function(freqs, freqs0=NULL) {
  if (!is.null(freqs0)) {
    L<-list(freqs, freqs0)
    freqs<-tapply(unlist(L), names(unlist(L)), sum)
  }
  return(table(freqs))
}

lognRF <- function(psi, n) {
  return(sum(log((rep(psi, n)+seq(0, n-1)))))
}
