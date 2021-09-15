#' Realized Covariance
#'
#' Function returns the vectorization of the lowest triangular of the Realized Covariance matrices for different frequencies.
#' @param data a \code{(T x N)} \code{xts} object containing the \code{N} price/return series over period \code{T}
#' @param freq a string defining the desired frequency for the Realized Covariance matrices between "daily", "monthly", "quarterly" or "yearly"
#' @param make.ret boolean, in case it is \code{TRUE} the data are converted in returns, \code{FALSE} otherwise
#' @param cholesky boolean, in case it is \code{TRUE} the Cholesky factors of the Realized Covariance matrices are calculated, \code{FALSE} by default
#' @return
#' \item{Realized Covariances}{a \eqn{M \times N(N+1)/2} matrix of realized covariances, where \emph{M} is the number of lower frequency data}
#' \item{Cholesky Factors (optional)}{a \eqn{M \times N(N+1)/2} matrix of Cholesky factors of the realized
#' covariance matrices, where \emph{M} is the number of lower frequency data}
#' \item{returns (optional)}{a \eqn{T \times N} matrix of returns, when \code{make.ret = TRUE}}
#' @references Andersen T.G., Bollerslev T., Diebold F.X. and Labys P. (2003), Modeling and Forecasting Realized Volatility. \emph{Econometrica}. 71: 579-625
#'
#' Barndorff-Nielsen O.E. and Shephard  N. (2002), Econometric analysis of realised volatility and its use in estimating stochastic volatility models \emph{Journal of the Royal Statistical Society}. 64(2): 253-280
#' @author Andrea Bucci
#' @importFrom methods is
#' @importFrom zoo as.yearqtr
#' @importFrom xts xts apply.daily apply.monthly apply.quarterly apply.yearly
#' @importFrom matrixcalc vech
#' @importFrom quantmod dailyReturn monthlyReturn yearlyReturn
#' @importFrom ks invvech
#' @importFrom zoo zoo
#' @importFrom lessR to
#' @export
#' @keywords RCOV
#' @examples
#' data(Sample5minutes)
#' rc <- rcov(Sample5minutes, freq = 'daily', cholesky = TRUE, make.ret = TRUE)
#' print(rc)
#'
rcov <- function(data, freq = c('daily', 'monthly', 'quarterly', 'yearly'), make.ret = TRUE, cholesky = FALSE){
freq <- match.arg(freq)
  if(is(data, 'xts') == FALSE){
    stop('Data should be of class "xts".')
  }
ncoly <- ncol(data)
nrowy <- nrow(data)
start_date <- index(data)[1]
end_date <- index(data)[nrowy]
if(freq == 'daily'){
  days <- as.Date(index(data), format = "%m/%d/%Y %I:%M:%S %p")
  days <- unique(days)
  nday <- length(days)
}

elapsed_months <- function(end_date, start_date){
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)+1
}
elapsed_quarters <- function(end_date, start_date){
  (as.yearqtr(end_date)-
     as.yearqtr(start_date))*4+1
}
elapsed_years <- function(end_date, start_date){
  as.numeric(format(end_date, format = '%Y'))-
    as.numeric(format(start_date, format = '%Y'))+1
}

makeReturns <- function (ts){
  l <- dim(ts)[1]
  col_names <- names(ts)
  x <- matrix(as.numeric(ts), nrow = l)
  x[(2:l), ] <- log(x[(2:l), ]) - log(x[(1:(l - 1)), ])
  x[1, ] <- rep(0, dim(ts)[2])
  x <- xts(x, order.by = index(ts))
  names(x) <- col_names
  return(x)
}
nmonth <- elapsed_months(end_date, start_date)
nquarter <- elapsed_quarters(end_date, start_date)
nyear <- elapsed_years(end_date, start_date)


retu <- matrix(ncol = ncoly, nrow = nrowy)
if(make.ret == TRUE){
    if(freq == 'daily'){
      for(j in 1:ncoly){
        retu[,j] <- makeReturns(data[,j])*100
      }
      crosspro <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
      for (i in 1:nrowy){
        cross <- retu[i,]%*%t(retu[i,])
        crosspro[i,] <- vech(cross)
      }
      cross1 <- zoo(crosspro, order.by = index(data))

      realized <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nday)
      for (k in 1:nday){
        realized[, k] <- apply.daily(cross1[,k], sum)
      }
      chol2 <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nday)
      for (j in 1:nday){
        rcovmat <- invvech(as.matrix(realized[j,]))
        chol1 <- t(chol(rcovmat, pivot = T))
        chol2[j,] <- vech(chol1)
      }
      realized <- as.xts(realized, order.by = days)
      chol2 <- as.xts(chol2, order.by = days)
    }
    if(freq == 'monthly'){
      for(j in 1:ncoly){
       retu[,j] <- dailyReturn(data[,j])*100
      }
      crosspro <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
      for (i in 1:nrowy){
        cross <- retu[i,]%*%t(retu[i,])
        crosspro[i,] <- vech(cross)
      }
      cross1 <- zoo(crosspro, order.by = index(data))

      realized <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nmonth)
      for (k in 1:(ncoly*(ncoly+1)/2)){
        realized[, k] <- apply.monthly(cross1[,k], sum)
      }
      chol2 <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nmonth)
      for (j in 1:nmonth){
        rcovmat <- invvech(as.matrix(realized[j,]))
        chol1 <- t(chol(rcovmat, pivot = T))
        chol2[j,] <- vech(chol1)
      }
      realized <- as.xts(realized, order.by = seq(start_date, end_date, by = 'month'))
      chol2 <- as.xts(chol2, order.by = seq(start_date, end_date, by = 'month'))
    }
    if(freq == 'quarterly'){
      for(j in 1:ncoly){
        retu[,j] <- dailyReturn(data[,j])*100
      }
      crosspro <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
      for (i in 1:nrowy){
        cross <- retu[i,]%*%t(retu[i,])
        crosspro[i,] <- vech(cross)
      }
      cross1 <- zoo(crosspro, order.by = index(data))

      realized <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nquarter)
      for (k in 1:(ncoly*(ncoly+1)/2)){
        realized[, k] <- apply.quarterly(cross1[,k], sum)
      }
      chol2 <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nquarter)
      for (j in 1:nquarter){
        rcovmat <- invvech(as.matrix(realized[j,]))
        chol1 <- t(chol(rcovmat, pivot = T))
        chol2[j,] <- vech(chol1)
      }
      realized <- as.xts(realized, order.by = seq(start_date, end_date, by = 'quarter'))
      chol2 <- as.xts(chol2, order.by = seq(start_date, end_date, by = 'quarter'))
    }
    if(freq == 'yearly'){
      for(j in 1:ncoly){
        retu[,j] <- dailyReturn(data[,j])*100
      }
      crosspro <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
      for (i in 1:nrowy){
        cross <- retu[i,]%*%t(retu[i,])
        crosspro[i,] <- vech(cross)
      }
      cross1 <- zoo(crosspro, order.by = index(data))

      realized <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nyear)
      for (k in 1:(ncoly*(ncoly+1)/2)){
        realized[, k] <- apply.yearly(cross1[,k], sum)
      }
      chol2 <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nyear)
      for (j in 1:nyear){
        rcovmat <- invvech(as.matrix(realized[j,]))
        chol1 <- t(chol(rcovmat, pivot = T))
        chol2[j,] <- vech(chol1)
      }
      realized <- as.xts(realized, order.by = seq(start_date, end_date, by = 'year'))
      chol2 <- as.xts(chol2, order.by = seq(start_date, end_date, by = 'year'))
    }
}else{
    retu <- data
    if(freq == 'daily'){
      crosspro <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
      for (i in 1:nrowy){
        cross <- retu[i,]%*%t(retu[i,])
        crosspro[i,] <- vech(cross)
      }
      cross1 <- zoo(crosspro, order.by = index(data))

      realized <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nday)
      for (k in 1:(ncoly*(ncoly+1)/2)){
        realized[, k] <- apply.daily(cross1[,k], sum)
      }
      chol2 <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nday)
      for (j in 1:nday){
        rcovmat <- invvech(as.matrix(realized[j,]))
        chol1 <- t(chol(rcovmat, pivot = T))
        chol2[j,] <- vech(chol1)
      }
      realized <- as.xts(realized, order.by = days)
      chol2 <- as.xts(chol2, order.by = days)
    }
    if(freq == 'monthly'){
      crosspro <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
      for (i in 1:nrowy){
        cross <- retu[i,]%*%t(retu[i,])
        crosspro[i,] <- vech(cross)
      }
      cross1 <- zoo(crosspro, order.by = index(data))

      realized <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nmonth)
      for (k in 1:(ncoly*(ncoly+1)/2)){
        realized[, k] <- apply.monthly(cross1[,k], sum)
      }
      chol2 <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nmonth)
      for (j in 1:nmonth){
        rcovmat <- invvech(as.matrix(realized[j,]))
        chol1 <- t(chol(rcovmat, pivot = T))
        chol2[j,] <- vech(chol1)
      }
      realized <- as.xts(realized, order.by = seq(start_date, end_date, by = 'month'))
      chol2 <- as.xts(chol2, order.by = seq(start_date, end_date, by = 'month'))
    }
    if(freq == 'quarterly'){
      crosspro <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
      for (i in 1:nrowy){
        cross <- retu[i,]%*%t(retu[i,])
        crosspro[i,] <- vech(cross)
      }
      cross1 <- zoo(crosspro, order.by = index(data))

      realized <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nquarter)
      for (k in 1:(ncoly*(ncoly+1)/2)){
        realized[, k] <- apply.quarterly(cross1[,k], sum)
      }
      chol2 <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nquarter)
      for (j in 1:nquarter){
        rcovmat <- invvech(as.matrix(realized[j,]))
        chol1 <- t(chol(rcovmat, pivot = T))
        chol2[j,] <- vech(chol1)
      }
      realized <- as.xts(realized, order.by = seq(start_date, end_date, by = 'quarter'))
      chol2 <- as.xts(chol2, order.by = seq(start_date, end_date, by = 'quarter'))
    }
    if(freq == 'yearly'){
      crosspro <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
      for (i in 1:nrowy){
        cross <- retu[i,]%*%t(retu[i,])
        crosspro[i,] <- vech(cross)
      }
      cross1 <- zoo(crosspro, order.by = index(data))

      realized <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nyear)
      for (k in 1:(ncoly*(ncoly+1)/2)){
        realized[, k] <- apply.yearly(cross1[,k], sum)
      }
      chol2 <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nyear)
      for (j in 1:nyear){
        rcovmat <- invvech(as.matrix(realized[j,]))
        chol1 <- t(chol(rcovmat, pivot = T))
        chol2[j,] <- vech(chol1)
      }
      realized <- as.xts(realized, order.by = seq(start_date, end_date, by = 'year'))
      chol2 <- as.xts(chol2, order.by = seq(start_date, end_date, by = 'year'))
    }}

names1 <- to('y', (ncoly*(ncoly+1)/2), same.size = FALSE)

colnames(realized) <- names1
colnames(chol2) <- names1

if(make.ret == T){
  if(cholesky == TRUE){
    results <- list(realized, chol2, retu)
    names(results) <- c('Realized Covariances', 'Cholesky Factors', 'returns')
  }else{
    results <- list(realized, retu)
    names(results) <- c('Realized Covariances', 'returns')
  }
}else{
  if(cholesky == TRUE){
  results <- list(realized, chol2)
  names(results) <- c('Realized Covariances', 'Cholesky Factors')
  }else{
  results <- list(realized)
  names(results) <- 'Realized Covariances'
  }
}
return(results)
}
