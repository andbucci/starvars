#' Print method for objects of class VLSTAR
#'
#' \sQuote{\code{print}} methods for class \sQuote{\code{VLSTAR}}.
#' @aliases print
#' @param x An object of class \sQuote{\code{VLSTAR}} obtained through \command{VLSTAR()}.
#' @param \dots further arguments to be passed to and from other methods
#' @references Terasvirta T. and Yang Y. (2014), Specification, Estimation and Evaluation of Vector Smooth Transition Autoregressive Models with Applications. \emph{CREATES Research Paper 2014-8}
#' @author Andrea Bucci
#' @return Print of VLSTAR results
#' @keywords VLSTAR
#' @seealso \code{\link{VLSTAR}}
#' @export
#'
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
