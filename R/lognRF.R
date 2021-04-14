#' Log of rising factorial psi(psi+1)...(psi+n-1)
#'
#' lognRF calculates log(psi)+log(psi+1)+...+log(psi+n-1), and it might be
#' good to include an optional floor for j, as in log(psi + j)
#' @export

lognRF <- function(psi, n) {
  return(sum(log((rep(psi, n)+seq(0, n-1)))))
}
