#' Fit the supervised classifier
#'
#' The below function trains the model according to training data x and labels y.
#' The output is a classwise list including the frequencies of the data, the MLE of
#' psi, the probability of the partition from the logarithm of the Ewens
#' sampling formula, and a coefficient that equals
#' $log(n!/(psi(psi+1)...(psi+n-1)).
#' Saving the frequencies, the training probability and the coefficient save us
#' from having to calculate these quantities again for the test data.
#' @param x data vector.
#' @param y label vector.
#' @return a  list of classwise lists, each with components
#' @return frequencies  Frequencies of values
#' @return psi          MLE of psi
#' @return tr.prob      Probability of the training data under partition exchangeability
#' @return coeffs       n!/(psi(psi+1)...(psi+n-1)
#'
#' @keywords Fit training data
#' @export
#' @examples
#' set.seed(111)
#' x1<-rPD(5000,10)
#' x2<-rPD(5000,100)
#' x<-c(x1,x2)
#' y1<-rep("1", 5000)
#' y2<-rep("2", 5000)
#' y<-c(y1,y2)
#' SPCuPE.fit(x,y)

SPCuPE.fit <- function(x, y) {
  #x<-Data
  #y<-labels
  y<-sapply(y, toString)
  results<-list()
  classes<-unique(y)


  for (yi in classes) {
    cdata<-x[which(yi==y)]
    cw.freqs<-table(cdata)
    cw.abunds<-abundances(cw.freqs)
    cw.PsiMLE<-MLEp(cw.abunds)$psi
    cw.probs<-dlPD(cw.abunds, cw.PsiMLE)
    cw.coeffs<-lfactorial(length(cdata))-lognRF(cw.PsiMLE, length(cdata))
    results[[yi]]<-list(frequencies=cw.freqs, psi=cw.PsiMLE, tr.prob=cw.probs, coeffs=cw.coeffs)
  }
  return(results)
}
