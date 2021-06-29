#' Fit the supervised classifier
#'
#' Trains the model according to training data x and labels y.
#' The output is a classwise list including the frequencies of the data, and the MLE of
#' psi.
#' @param x data vector, or matrix with rows as data points and columns as features.
#' @param y training label vector.
#' @return If x is multidimensional, each list described below is returned for each dimension.
#' @return Returns a list of classwise lists, each with components:
#' @return frequencies: the frequencies of values in the class.
#' @return psi: the estimate of psi for the class.
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
#' fit<-SPEC.fit(x,y)
#'
#' ##With multidimensional x:
#' set.seed(111)
#' x1<-cbind(rPD(5000,10),rPD(5000,50))
#' x2<-cbind(rPD(5000,100),rPD(5000,500))
#' x<-rbind(x1,x2)
#' y1<-rep("1", 5000)
#' y2<-rep("2", 5000)
#' y<-c(y1,y2)
#' fit<-SPEC.fit(x,y)

SPEC.fit <- function(x, y) {
  if(length(dim(x))>1){
    feature_wise_results<-apply(x, 2, function(z) SPEC.fit(z, y))
    return(feature_wise_results)
  }

  y<-sapply(y, toString)
  results<-list()
  classes<-unique(y)


  for (yi in classes) {
    cdata<-x[which(yi==y)]
    cw.freqs<-table(cdata)
    cw.abunds<-abundances(cw.freqs)
    cw.PsiMLE<-MLEp(cw.abunds)$psi
    results[[yi]]<-list(frequencies=cw.freqs, psi=cw.PsiMLE)
  }
  return(results)
}





#' Marginally predicted labels of the test data given training data classification.
#'
#' tMarLab classifies the test data x based on the training data object.
#' The test data is considered i.i.d. and to have arrived sequentially. Thus, each
#' data point is classified one by one.
#' @param training A training data object from the function SPEC.fit.
#' @param x Test data vector or matrix with rows as data points and columns as features.
#' @return A vector of predicted labels for test data x.
#' @keywords Marginal classifier
#' @export
#' @references #' The classification algorithm is adapted from
#' [Corander, J., Cui, Y., Koski, T., and Siren, J.: Have I seen you before? Principles of Bayesian predictive classification revisited. Springer, Stat. Comput. 23, (2011), 59–73.]
#'  (https://doi.org/10.1007/s11222-011-9291-7)
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
#' fit<-SPEC.fit(x,y)
#'
#' tM<-tMarLab(fit, t)
#'
#' ##With multidimensional x:
#' set.seed(111)
#' x1<-cbind(rPD(5500,10),rPD(5500,50))
#' x2<-cbind(rPD(5500,100),rPD(5500,500))
#' test.ind1<-sample.int(5500,500)
#' test.ind2<-sample.int(5500,500)
#' x<-rbind(x1[-test.ind1,],x2[-test.ind2,])
#' y1<-rep("1", 5000)
#' y2<-rep("2", 5000)
#' y<-c(y1,y2)
#' fit<-SPEC.fit(x,y)
#' t1<-x1[test.ind1,]
#' t2<-x2[test.ind2,]
#' t<-rbind(t1,t2)
#'
#' tM<-tMarLab(fit, t)



tMarLab <- function(training, x) {
  pred<-c()
  if(!is.null(training[[1]]$psi)){training<-list(training)}
  nfeatures<-length(training)
  classes<-names(training[[1]])
  x<-matrix(x, ncol = nfeatures)
  for (i in 1:length(x[,1])) { ####problem: how to for entire row
    probs<-cbind(classes, rep(1.0, length(classes))) #vector to collect prob of each class
    j<-0 #class index for prob vector
    xi<-x[i,]
    for (class in classes) {
      j<-j+1
      for(feat in 1:nfeatures) {
          xd<-toString(xi[feat])
        if (xd %in% names(training[[feat]][[class]][["frequencies"]])) {
          probs[j,2]<- as.double(probs[j,2]) * training[[feat]][[class]][["frequencies"]][[xd]] /
            (sum(training[[feat]][[class]][["frequencies"]]) + training[[feat]][[class]][["psi"]])
        } else {
          probs[j,2]<- as.double(probs[j,2]) * training[[feat]][[class]][["psi"]] /
            (sum(training[[feat]][[class]][["frequencies"]]) + training[[feat]][[class]][["psi"]])
        }
      }
    }
    pred<-c(pred, probs[which.max(probs[,2]),1], use.names=FALSE)

  }
  return(pred)
}





#' Simultaneously predicted labels of the test data given the training data classification.
#'
#' The simultaneous case:
#' The test data are first labeled with the marginal classifier. The simultaneous
#' classifier then iterates over all test data, assigning each a label by finding
#' the maximum predictive probability given the current classification structure of
#' the test data as a whole. This is repeated until the classification structure
#' doesn't change after iterating over all data.
#' The classification algorithm is adapted from
#' [Corander, J., Cui, Y., Koski, T., and Siren, J.: Have I seen you before? Principles of Bayesian predictive classification revisited. Springer, Stat. Comput. 23, (2011), 59–73.] (https://doi.org/10.1007/s11222-011-9291-7)
#' @param training A training data object from the function SPEC.fit.
#' @param x Test data vector or matrix with rows as data points and columns as features.
#' @return A vector of predicted labels for test data x.
#' @keywords Simultaneous classifier
#' @usage tSimLab(training, x)
#' @export
#' @references The classification algorithm is adapted from
#' [Corander, J., Cui, Y., Koski, T., and Siren, J.: Have I seen you before? Principles of Bayesian predictive classification revisited. Springer, Stat. Comput. 23, (2011), 59–73.]
#'  (https://doi.org/10.1007/s11222-011-9291-7)
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
#' fit<-SPEC.fit(x,y)
#'
#' tS<-tSimLab(fit, t)
#'
#' ##With multidimensional x:
#' set.seed(111)
#' x1<-cbind(rPD(5500,10),rPD(5500,50))
#' x2<-cbind(rPD(5500,100),rPD(5500,500))
#' test.ind1<-sample.int(5500,500)
#' test.ind2<-sample.int(5500,500)
#' x<-rbind(x1[-test.ind1,],x2[-test.ind2,])
#' y1<-rep("1", 5000)
#' y2<-rep("2", 5000)
#' y<-c(y1,y2)
#' fit<-SPEC.fit(x,y)
#' t1<-x1[test.ind1,]
#' t2<-x2[test.ind2,]
#' t<-rbind(t1,t2)
#'
#' tS<-tSimLab(fit, t)



tSimLab <- function(training, x) {
  if(!is.null(training[[1]]$psi)){training<-list(training)}
  classes<-names(training[[1]])
  nfeatures<-length(training)
  x<-matrix(x, ncol = nfeatures)

  #pred<-cbind(x,tMarLab(training,x)) #using the new tMarLab. pred is the iterator
  #predS <- cbind(rep(NULL, length(x)), rep(NULL, length(x))) # predS is the previous value of pred

  pred<-tMarLab(training, x)
  predS<-rep(NULL, length(pred))

  #while (sum(pred[,2]==predS[,2])!=length(x[,1])) {
  while (sum(pred==predS)!=length(x[,1])) {
    predS<-pred
    for (i in 1:length(x[,1])) {
      probsByClass<-c()
      for (class in classes) {
        pred[i] <- NaN
        prob<-1
        #combined frequencies of training and test data in the same class
        #add a line that takes the currently predicted test item out of test data
        for (feat in 1:nfeatures) {

          L<-list(table(x[which(pred==class),feat]), training[[feat]][[class]][["frequencies"]])#combined frequencies of training and test data in the same class
          freqs<-tapply(unlist(L), names(unlist(L)), sum)

          xid<-toString(x[i,feat])
          if (xid %in% names(freqs)) {
            prob<-prob*freqs[[xid]] /
              (sum(freqs) + training[[feat]][[class]][["psi"]])
          } else {
            prob<-prob*training[[feat]][[class]][["psi"]] /
              (sum(freqs) + training[[feat]][[class]][["psi"]])
          }
        }

        probsByClass <- append(probsByClass, prob)
      }
      pred[i] <- classes[which.max(probsByClass)]
    }
  }

  return(pred)
}


