#' Paired test based on psi=0.
#'


PE.prdtest<-function(x1,x2) {
  x=abs(x1-x2)
  mle=MLEp(table(table(x)))

}

MLEp.0ci<- function(n) {


  #n<-sum(as.integer(names(abund))*abund)
  target<- exp(log(0.05) - lfactorial(n-1))
  psi<-0.1
  asum<- exp(log(psi) - lognRF(psi, n))
  last<-psi*2
  print(asum-target)
  while(abs(asum-target)>.Machine$double.xmin) {


    if (asum>target && last==psi*2) {
      last<-psi
      psi<-2*psi
    } else if (asum>target){
      psi1<-psi
      psi<- psi + abs(psi-last)/2
      last<-psi1
    } else if (asum<target && last==psi*2){
      last<-psi
      psi<-0.5*psi
    } else {
      psi1<-psi
      psi<- psi - abs(psi-last)/2
      last<-psi1
    }

    if(last==psi) {break}

    asum<- exp(log(psi) - lognRF(psi, n))   #sum(psi/(psi+seq(0,n-1)))
  }

  return(psi)
}
