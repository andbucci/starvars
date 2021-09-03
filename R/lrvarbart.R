#' Long-run variance using Bartlett kernel
#'
#' Function returns the long-run variance of a time series, relying on the Bartlett kernel.
#' The window size of the kernel is the cube root of the sample size.
#' @param x a \code{(T x 1)} vector containing the time series over period \code{T}
#' @return
#' \item{lrv}{long-run variance}
#' \item{return}{bandwidth size of the window}
#' @references Hamilton J. D. (1994), Time Series Analysis. \emph{Princeton University Press}; Tsay R.S. (2005), Analysis of Financial Time Series. \emph{John Wiley & SONS}
#' @author Andrea Bucci
#' @keywords VLSTAR
#' @export
#' @importFrom stats acf
#' @examples
#' data(Realized)
#' lrvarbart(Realized[,1])
#'
lrvarbart <- function(x){
  x <- as.xts(x)
    N <- length(x)
    cc = 1.4
    K = N^(1/3)
    ac <- acf(x, plot = FALSE, lag.max = floor(2 * N^(2/3)))$acf
    vc <- cc * sqrt(log(N, base = 10)/N)
    for (i in 1:floor(N^(2/3))){
      if(max(abs(ac[i + (1:K)])) < vc){
        break
      }
    }
    w <- ((2 * i):1)/(2 * i)
    ac <- acf(x, plot = FALSE, type = "covariance", lag.max = 2*i)$acf[1:(2 * i + 1)]
    asy <- ac[1] + 2 * sum(ac[2:(2 * i + 1)] * w)
    erg <- list(lrv = asy, bandwidth = i)
    return(erg)
}

