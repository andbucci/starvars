print.VLSTAR <- function(object, ...) {
  digits = 3
  cat(paste("\nVLSTAR model Estimation through", object$method, "\n"))
  order = (object$m-1):object$m
  ord.coef = list()
  for(i in 1:object$m){
    ord.coef[[i]] = as.data.frame(object$Bhat[grep(paste("m_ ", order[i], sep=''), rownames(object$Bhat)),])
  }
  names(ord.coef) = paste("Regime ", order, sep='')
  gammaCoef <- object$Gammac[,1]
  cCoef <- object$Gammac[,2]
  
  cat("Coefficients:\n")
  print(ord.coef)
  cat("\nSmoothing parameter: gamma =", format(gammaCoef, digits=4),"\n")
  cat("\nThreshold")
  cat("\nValue:", format(cCoef, digits=4), "\n")
  invisible(object)
  #NextMethod('print')
}