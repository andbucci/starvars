
#' @S3method print VLSTARnls
print.VLSTAR <- function(x, ...) {
  #NextMethod(...)
  cat("\nVLSTAR model Estimation through Nonlinear Least Squares\n")
  order.L <- (x$m-1)
  order.H <- x$m
  lowCoef <- x$Bhat[grep(paste("m_ ", order.L, sep=''), rownames(x$Bhat))]
  highCoef <- x$Bhat[grep(paste("m_ ", order.H, sep=''), rownames(x$Bhat))]
  gammaCoef <- x$Cgamma[,1]
  cCoef <- x$Cgamma[,2]
  
  cat("Coefficients:\n")
  cat("Low regime:\n")
  print(lowCoef, ...)
  cat("\nHigh regime:\n")
  print(highCoef, ...)
  cat("\nSmoothing parameter: gamma =", format(gammaCoef, digits=4),"\n")
  cat("\nThreshold")
  cat("\nValue:", format(cCoef, digits=4), "\n")
  invisible(x)
}