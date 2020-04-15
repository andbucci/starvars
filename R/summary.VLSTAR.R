#' @S3method summary VLSTAR
summary.VLSTAR<-function(object,...){
  #NextMethod(...)
  x<-object
  k<-ncol(x$Data[[2]])
  t<-nrow(x$Data[[1]])
  p<-x$p
  x$T <- nrow(x$yoriginal)
  Z<-t(as.matrix(tail.matrix(x$Data[[1]])))
  x$npar <- k*ncol(x$Data[[1]])*x$m + 2*(x$m-1)*ncol(x$Data[[1]])
  
  ## export results
  x$coefficients<-as.list(as.data.frame(x$Bhat))
  x$StDev<-as.list(as.data.frame(x$StDev))
  x$Pvalues<-as.list(as.data.frame(x$pval))
  x$Tvalues<-as.list(as.data.frame(x$ttest))
  ab<-list()
  symp<-list()
  stars<-list()
  for(i in 1:length(x$Pvalues)){
    symp[[i]] <- symnum(x$Pvalues[[i]], corr=FALSE,cutpoints = c(0,  .001,.01,.05, .1, 1), symbols = c("***","**","*","."," "))
    stars[[i]]<-matrix(symp[[i]], nrow=length(x$Pvalues[[i]]))
    ab[[i]]<-matrix(paste(x$coefficients[[i]],"(", x$StDev[[i]],")",stars[[i]], sep=""), nrow=length(x$StDev[[i]]))
    dimnames(ab[[i]])<-dimnames(x$coefficients[[1]])
  }
  attributes(ab)<-attributes(x$coefficients)
  x$stars<-stars
  x$starslegend<-symp[[1]]
  x$bigcoefficients<-ab
  x$aic<-x$AIC
  x$bic<-x$BIC
  x$t <- t
  x$k <- k
  class(x) <- 'summary.VLSTAR'
  return(x)
}