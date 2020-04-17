lrvarbart <- function(x){
  x <- xts::as.xts(x)
    N <- length(x)
    cc = 1.4
    K = N^(1/3)
    ac <- stats::acf(x, plot = F, lag.max = floor(2 * N^(2/3)))$acf
    vc <- cc * sqrt(log(N, base = 10)/N)
    for (i in 1:floor(N^(2/3))){
      if(max(abs(ac[i + (1:K)])) < vc){
        break
      }
    }
    w <- ((2 * i):1)/(2 * i)
    ac <- stats::acf(x, plot = F, type = "covariance", lag.max = 2*i)$acf[1:(2 * i + 1)]
    asy <- ac[1] + 2 * sum(ac[2:(2 * i + 1)] * w)
    erg <- list(lrv = asy, bandwidth = i)
    return(erg)
}

