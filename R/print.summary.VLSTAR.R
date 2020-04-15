#' @S3method print summary.VLSTAR
print.summary.VLSTAR<-function(x,digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"),...){
  coeftoprint<-list()
  for(i in 1:length(x$bigcoefficients)){
    a<-myformat(x$coefficients[[i]], digits)
    b<-myformat(x$StDev[[i]], digits)
    aic1 <- round(x$AIC,2)
    bic1 <- round(x$BIC,2)
    if(getOption("show.signif.stars"))
      stars<-x$stars[[i]]
    else
      stars<-NULL
    coeftoprint[[i]]<-matrix(c(paste(a,"(", b,")",stars, sep=""), '', round(x$Cgamma[i,1],4), round(x$Cgamma[i,2],4),'',
                               aic1[i], bic1[i]), nrow=(length(x$StDev[[1]])+6))
    rownames(coeftoprint[[i]])<- c(rownames(x$Bhat),"---",'gamma', 'c',"---", 'AIC', 'BIC')
  }
  names(coeftoprint) <- colnames(x$Bhat)
  cat("Model VLSTAR with ", x$m, " regimes\n", sep ='')
  cat("\nFull sample size:",x$T, "\tEnd sample size:", x$t)
  cat("\nNumber of variables:", x$k,"\tNumber of estimated parameters:", x$npar)
  #cat("\nAIC",x$aic)
  #cat("\nBIC", x$bic)
  cat("\nMultivariate log-likelihood:", x$MultiLL,"\n\n")
  print(noquote(coeftoprint))
  if (signif.stars)
    cat("---\nSignif. codes: ", attr(x$starslegend, "legend"), "\n")
  #cat("\nThreshold values:",x$Cgamma[,2])
}