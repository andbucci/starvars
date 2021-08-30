print.VLSTAR <- function(x, ...) {
  digits = 3
  cat(paste("\nVLSTAR model Estimation through", x$method, "\n"))
  order = (x$m-1):x$m
  ord.coef = list()
  for(i in 1:x$m){
    ord.coef[[i]] = as.data.frame(x$Bhat[grep(paste("m_ ", order[i], sep=''), rownames(x$Bhat)),])
  }
  names(ord.coef) = paste("Regime ", order, sep='')
  gammaCoef <- x$Gammac[,1]
  cCoef <- x$Gammac[,2]

  cat("Coefficients:\n")
  print(ord.coef)
  cat("\nSmoothing parameter: gamma =", format(gammaCoef, digits=4),"\n")
  cat("\nThreshold")
  cat("\nValue:", format(cCoef, digits=4), "\n")
  invisible(x)
  #NextMethod('print')
}
