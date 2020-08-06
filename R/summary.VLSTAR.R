#' @S3method summary VLSTAR
summary.VLSTAR<-function(object,...){
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
  NextMethod('print.summary')
}

#' @S3method print summary.VLSTAR
print.summary.VLSTAR<-function(x,...){
  digits = 3
  coeftoprint<-list()
  myformat<-function(x,digits, toLatex=FALSE){
    r<-x
    littlex<-abs(x)<10^-(digits)
    r[!littlex]<-formatC(x[!littlex],digits=digits, format="f")
    r[littlex]<-format(x[littlex],digits=min(digits,2), scientific=TRUE)
    if(toLatex)
      r<-gsub("(e.*)","slashtext{\\1}",r)
    if(class(x)=="numeric")
      return(noquote(r))
    if(class(x)=="matrix")
      return(matrix(noquote(r), ncol=ncol(x), nrow=nrow(x)))
  }
  for(i in 1:length(x$bigcoefficients)){
    a<-myformat(x$coefficients[[i]], digits)
    b<-myformat(x$StDev[[i]], digits)
    c<-myformat(x$Pvalues[[i]], digits)
    aic1 <- round(x$AIC,2)
    bic1 <- round(x$BIC,2)
    stars1<-x$stars[[i]]
    coeftop <- cbind(a, b, c, stars1)
    colnames(coeftop) <- c('Estimate', 'Std. Error', 'p-value', '')
    coeftoprint[[i]] <- rbind(coeftop, rep('', 4), c(round(x$Gammac[i,1],4),'','',''), c(round(x$Gammac[i,2],4),'','',''),
                              rep('', 4), c(aic1[i],'','',''), c(bic1[i], '','',''))
    rownames(coeftoprint[[i]])<- c(rownames(x$Bhat),"---",'gamma', 'c',"---", 'AIC', 'BIC')
  }
  names(coeftoprint) <- colnames(x$Bhat)
  cat("Model VLSTAR with ", x$m, " regimes\n", sep ='')
  cat("\nFull sample size:",x$T)
  cat("\nNumber of variables used as covariates:", x$k,"\tNumber of estimated parameters:", x$npar)
  #cat("\nAIC",x$aic)
  #cat("\nBIC", x$bic)
  cat("\nMultivariate log-likelihood:", x$MultiLL,"\n\n")
  cat('\nCoefficients:')
  print(noquote(coeftoprint))
  cat("=================================\n")
  cat("\n Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
}
