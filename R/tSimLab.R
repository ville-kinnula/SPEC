#' Simultaneously predicted labels of the test data given the training data classification.
#'
#' The simultaneous case:
#' The test data are first labeled with the marginal classifier. The simultaneous
#' classifier then iterates over all test data, assigning each a label by finding
#' the maximum predicting probability given the current classification structure of
#' the test data as a whole. This is repeated until the classification structure
#' doesn't change after iterating over all data.
#' The classification algorithm is from
#' [Corander, J., Cui, Y., Koski, T., and Siren, J.: Have I seen you before? Principles of Bayesian predictive classification revisited. Springer, Stat. Comput. 23, (2011), 59–73.] (https://doi.org/10.1007/s11222-011-9291-7)
#' @param training A training data object from the function SPCuPE.fit.
#' @param x Test data vector.
#' @return A vector of predicted labels for test data x.
#' @keywords Marginal classifier
#' @usage tSimLab(training, x)
#' @export
#' @examples
#' set.seed(111)
#' x1<-rPD(11500,10)
#' x2<-rPD(11500,1000)
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
#' tS<-tSimLab(fit, t)

tSimLab <- function(training, x) {
  classes<-names(training)
  pred<-cbind(x,tMarLab(training,x))
  predS <- cbind(rep(NULL, length(x)), rep(NULL, length(x)))
  while (sum(pred[,2]==predS[,2])!=length(x)) {
    predS<-pred
    for (i in 1:length(x)) {
      probsByClass<-c()
      for (class in classes) {
        pred[i,2] <- class
        probsByClass <- append(probsByClass,
                               sum(
                                 sapply(
                                   classes,
                                   function (y) dPD.pred(
                                     pred[which(pred[,2]==y),1],
                                     training[[y]][["psi"]],
                                     training[[y]][["frequencies"]],
                                     training[[y]][["coeffs"]],
                                     training[[y]][["tr.prob"]])
                                 )))
      }
      pred[i,2] <- classes[which.max(probsByClass)]
    }
  }

  return(pred[,2])
}


#' Simultaneously predicted labels of the test data given the training data classification.
#'
#' The simultaneous case:
#' The test data are first labeled with the marginal classifier. The simultaneous
#' classifier then iterates over all test data, assigning each a label by finding
#' the maximum predicting probability given the current classification structure of
#' the test data as a whole. This is repeated until the classification structure
#' doesn't change after iterating over all data.
#' The classification algorithm is from
#' [Corander, J., Cui, Y., Koski, T., and Siren, J.: Have I seen you before? Principles of Bayesian predictive classification revisited. Springer, Stat. Comput. 23, (2011), 59–73.] (https://doi.org/10.1007/s11222-011-9291-7)
#' @param training A training data object from the function SPCuPE.fit.
#' @param x Test data vector.
#' @return A vector of predicted labels for test data x.
#' @keywords Marginal classifier
#' @usage tSimLab(training, x)
#' @export
#' @examples
#' set.seed(111)
#' x1<-rPD(11500,10)
#' x2<-rPD(11500,1000)
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
#' tS<-tSimLab2(fit, t)


tSimLab2 <- function(training, x) {
  classes<-names(training)
  pred<-cbind(x,tMarLab2(training,x)) #using the new tMarLab. pred is the iterator
  predS <- cbind(rep(NULL, length(x)), rep(NULL, length(x))) # predS is the previous value of pred
  while (sum(pred[,2]==predS[,2])!=length(x)) {
    predS<-pred
    for (i in 1:length(x)) {
      probsByClass<-c()
      for (class in classes) {
        pred[i,2] <- class
        #combined frequencies of training and test data in the same class
        L<-list(table(pred[which(pred[,2]==class),1]), training[[class]][["frequencies"]])#combined frequencies of training and test data in the same class
        freqs<-tapply(unlist(L), names(unlist(L)), sum)

        prob<-0
        if (toString(x[i]) %in% names(freqs)) { #BUG? FIXED?
          prob<-freqs[[toString(x[i])]] /
            (sum(freqs) + training[[class]][["psi"]])
        } else {
          prob<-training[[class]][["psi"]] /
            (sum(freqs) + training[[class]][["psi"]])
        }

        probsByClass <- append(probsByClass, prob)
      }
      pred[i,2] <- classes[which.max(probsByClass)]
    }
  }

  return(pred[,2])
}
