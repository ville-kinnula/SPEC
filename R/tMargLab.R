#' Marginally predicted labels of the test data given training data classification.
#'
#' tMarLab classifies the test data x based on the training data object.
#' The test data is considered i.i.d. and to have arrived sequentially. Thus, each
#' data point is classified one by one.
#' The classification algorithm is from
#' [Corander, J., Cui, Y., Koski, T., and Siren, J.: Have I seen you before? Principles of Bayesian predictive classification revisited. Springer, Stat. Comput. 23, (2011), 59–73.]
#'  (https://doi.org/10.1007/s11222-011-9291-7)
#' @param training A training data object from the function SPCuPE.fit.
#' @param x Test data vector.
#' @return A vector of predicted labels for test data x.
#' @keywords Marginal classifier
#' @export
#' @examples
#' set.seed(111)
#' x1<-rPD(10500,10)
#' x2<-rPD(10500,1000)
#' test.ind1<-sample.int(10500,500)
#' test.ind2<-sample.int(10500,500)
#' x<-c(x1[-test.ind1],x2[-test.ind2])
#' y1<-rep("1", 10000)
#' y2<-rep("2", 10000)
#' y<-c(y1,y2)
#'
#' t1<-x1[test.ind1]
#' t2<-x2[test.ind2]
#' t<-c(t1,t2)
#'
#' fit<-SPCuPE.fit(x,y)
#'
#' tM<-tMarLab(fit, t)
#'
#'
#' v<-c(x1[1:10000],x2[1:10000])
#' z<-c(x1[10001:10500],x2[10001:10500])
#'
#' fit2<-SPCuPE.fit(v,y)
#' tM2<-tMarLab(fit2,z)


tMarLab <- function(training, x) {
  pred <- c()
  classes<-names(training)

  for (i in 1:length(x)) {
    probs<-cbind(classes, rep(0, length(classes)))
    j<-0
    for (class in classes) {
      j<-j+1
      probs[j,2]<- dPD.pred(x[i],
                            training[[class]][["psi"]],
                            training[[class]][["frequencies"]],
                            training[[class]][["coeffs"]],
                            training[[class]][["tr.prob"]])
    }
    pred<-c(pred, probs[which.max(probs[,2]),1], use.names=FALSE)
  }
  return(pred)
}
