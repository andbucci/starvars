coef.VLSTAR <-
  function(object, ...){
    ## export results
    coefs<-as.list(as.data.frame(object$Bhat))
    gamma = object$Gammac[,1]
    c = object$Gammac[,2]
    names(coefs) <- colnames(object$yoriginal)
    coeftoprint = list(coefs, gamma, c)
    names(coeftoprint) = c('coefficients', 'gamma', 'c')
    return(coefs)
  }
