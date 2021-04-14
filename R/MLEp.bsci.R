#' Bootstrap confidence interval for the MLE of psi
#'
#' A bootstrapped confidence interval for the Maximum Likelihood Estimate for psi.
#' @param x A data vector.
#' @export
#' @return Lower and upper bounds of the confidence interval.


MLEp.bsci<-function(x) {
  n_bootstrap<-1000
  n<-length(x)
  bootsrap_psis<-c()
  for (i in 1:n_bootstrap) {
    bootsrap_psis<-append(bootsrap_psis, MLEp(table(table(sample(x,n,replace=TRUE))))$psi)
  }
  return(quantile(bootsrap_psis,c(0.025,0.975)))
}
