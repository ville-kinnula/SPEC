#' Predictive probability for a marginal test data classification.
#'
#' @param x.tr Training data vector.
#' @param y.tr Training data labels.
#' @param x.te Test data vector.
#' @param y.te Predicted test data labels
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
#' mar.prob(x,y,t,tM)


mar.prob <- function(x.tr, y.tr, x.te, y.te) {
  y.tr<-sapply(y.tr, toString)
  y.te<-sapply(y.te, toString)
  training<-SPCuPE.fit(x.tr,y.tr)
  te.data<-cbind(x.te,y.te)

  prob<-exp(sum(apply(te.data, 1, function(x) dPD.pred(x[1],
                                                       training[[x[2]]][["psi"]],
                                                       training[[x[2]]][["frequencies"]],
                                                       training[[x[2]]][["coeffs"]],
                                                       training[[x[2]]][["tr.prob"]]))))
  return(prob)
}
