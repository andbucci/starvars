coef.VLSTAR <-
  function(object, ...){
    x<-object
    ## export results
    x$coefficients<-as.list(as.data.frame(x$Bhat))
    x$StDev<-as.list(as.data.frame(x$StDev))
    x$Pvalues<-as.list(as.data.frame(x$pval))
    x$Tvalues<-as.list(as.data.frame(x$ttest))
    digits = 3
    coeftoprint<-list()
    for(i in 1:length(x$coefficients)){
      a<-signif(x$coefficients[[i]], 3)
      b<-signif(x$StDev[[i]], 3)
      c<-signif(x$Pvalues[[i]], 3)
      coeftop <- as.data.frame(cbind(a, b, c))
      colnames(coeftop) <- c('Estimate', 'Std. Error', 'p-value')
      if(x$singlecgamma == TRUE){
        coeftoprint[[i]] <- rbind(coeftop, rep('', 3), c(round(x$Cgamma[1,1],3),'',''), c(round(x$Cgamma[1,2],3),'',''))
      } else{
        coeftoprint[[i]] <- rbind(coeftop, rep('', 3), c(round(x$Gammac[i,1],3),'',''), c(round(x$Gammac[i,2],3),'',''))
        }
      rownames(coeftoprint[[i]])<- c(rownames(x$Bhat),"---",'gamma', 'c')
    }
    names(coeftoprint) <- colnames(object$yoriginal)
    return(coeftoprint)
  }
