% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Parameter_estimation_and_hypothesis_testing.R
\name{MLEp}
\alias{MLEp}
\title{Maximum Likelihood Estimate for psi}
\usage{
MLEp(abund)
}
\arguments{
\item{abund}{An abundance vector.}
}
\value{
Returns a list that contains the estimate "psi" and "Asymptotic confidence interval".
The confidence interval is based on the asymptotic distribution of maximum likelihood estimators.
}
\description{
Numerically searches for the MLE of psi given the frequencies of the
frequencies of a data vector, called an abundance vector.
}
\details{
Numerically searches for the MLE of psi as the root of equation
\deqn{K=\sum_{i=1}^n\psi/(\psi+i-1),} where \eqn{K} is the observed number of
different species in the sample. The right side of the equation is strictly
increasing when \eqn{\psi>0}, so a binary search is used to find the root.
An accepted \eqn{\psi} sets value of the right side
of the equation within R's smallest possible value of the actual value of \eqn{K}.
}
\examples{
##Find the MLE of psi of the vector (1,2,2).
##The frequencies of the frequencies of the data vector are given as input:
MLEp(table(table(c(1,2,2))))

##Find the MLE of psi of a sample from the Poisson-Dirichlet distribution:
set.seed(1000)
x<-rPD(n=10000, psi=100)
MLEp(table(table(x)))
}
\references{
W.J. Ewens, The sampling theory of selectively neutral alleles, Theoretical Population Biology, Volume 3, Issue 1,
1972, Pages 87-112, ISSN 0040-5809, \url{https://doi.org/10.1016/0040-5809(72)90035-4}.
}
\keyword{estimate}
\keyword{likelihood}
\keyword{maximum}
\keyword{psi}
