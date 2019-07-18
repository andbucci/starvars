rcov <- function(data, freq, ret = FALSE, cholesky = FALSE){
require(quantmod)
require(highfrequency)
require(ks)
require(lessR)
require(lubridate)
require(xts)
  if(!freq %in% c('daily', 'monthly', 'quarterly', 'yearly')){
    stop('Please, provide a correct frequency for the realized covariance matrix!')
  }
  if(is(data, 'xts') == F){
    stop('Data should be of class "xts".')
  }
ncoly <- ncol(data)
nrowy <- nrow(data)
start_date <- index(data)[1]
end_date <- index(data)[nrowy]
freq <- match.arg(freq)
if(freq == 'daily'){
  days <- as.Date(index(data), format = "%m/%d/%Y %I:%M:%S %p")
  days <- unique(days)
}

elapsed_months <- function(end_date, start_date){
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)+1
}
elapsed_quarters <- function(end_date, start_date){
  (as.yearqtr(end_date)-
     as.yearqtr(start_date))*4
}
elapsed_years <- function(end_date, start_date){
  year(end_date)-
    year(start_date)
}
nday <- length(days)
nmonth <- elapsed_months(end_date, start_date)
nquarter <- elapsed_quarters(end_date, start_date)
nyear <- elapsed_years(end_date, start_date)


retu <- matrix(ncol = ncoly, nrow = nrowy)
  if(ret == TRUE){
    if(freq == 'daily'){
      realized <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nday)
      for (k in 1:nday){
        realized[k,] <- vech(rCov(rdata = data[as.character(days[k])], makeReturns = T))
      }
      chol2 <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nday)
      for (j in 1:nday){
        rcovmat <- invvech(as.matrix(realized[j,]))
        chol1 <- t(chol(rcovmat))
        chol2[j,] <- vech(chol1)
      }
    }
    if(freq == 'monthly'){
      for(j in 1:ncoly){
       retu[,j] <- dailyReturn(data[,j])*100
      }
      crosspro <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
      for (i in 1:nrowy){
        cross <- t(retu[i,])%*%retu[i,]
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
        chol1 <- t(chol(rcovmat))
        chol2[j,] <- vech(chol1)
      }
    }
    if(freq == 'quarterly'){
      for(j in 1:ncoly){
        retu[,j] <- dailyReturn(data[,j])*100
      }
      crosspro <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
      for (i in 1:nrowy){
        cross <- t(retu[i,])%*%retu[i,]
        crosspro[i,] <- vech(cross)
      }
      cross1 <- zoo(crosspro, order.by = index(data))

      realized <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nquarter)
      for (k in 1:(ncoly*(ncoly+1)/2)){
        realized[, k] <- apply.quarterly(cross1[,k], sum)
      }
      chol2 <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nquarter)
      for (j in 1:nmonth){
        rcovmat <- invvech(as.matrix(realized[j,]))
        chol1 <- t(chol(rcovmat))
        chol2[j,] <- vech(chol1)
      }
    }
    if(freq == 'yearly'){
      for(j in 1:ncoly){
        retu[,j] <- monthlyReturn(data[,j])*100
      }
      crosspro <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
      for (i in 1:nrowy){
        cross <- t(retu[i,])%*%retu[i,]
        crosspro[i,] <- vech(cross)
      }
      cross1 <- zoo(crosspro, order.by = index(data))

      realized <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nyear)
      for (k in 1:(ncoly*(ncoly+1)/2)){
        realized[, k] <- apply.quarterly(cross1[,k], sum)
      }
      chol2 <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nyear)
      for (j in 1:nyear){
        rcovmat <- invvech(as.matrix(realized[j,]))
        chol1 <- t(chol(rcovmat))
        chol2[j,] <- vech(chol1)
      }}}else{
    retu <- data
    if(freq == 'daily'){
      realized <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nday)
      for (k in 1:nday){
        realized[k,] <- vech(rCov(rdata = data[as.character(days[k])], makeReturns = F))
      }
      chol2 <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nday)
      for (j in 1:nday){
        rcovmat <- invvech(as.matrix(realized[j,]))
        chol1 <- t(chol(rcovmat))
        chol2[j,] <- vech(chol1)
      }
    }
    if(freq == 'monthly'){
      crosspro <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
      for (i in 1:nrowy){
        cross <- t(retu[i,])%*%retu[i,]
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
        chol1 <- t(chol(rcovmat))
        chol2[j,] <- vech(chol1)
      }
    }
    if(freq == 'quarterly'){
      crosspro <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
      for (i in 1:nrowy){
        cross <- t(retu[i,])%*%retu[i,]
        crosspro[i,] <- vech(cross)
      }
      cross1 <- zoo(crosspro, order.by = index(data))

      realized <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nquarter)
      for (k in 1:(ncoly*(ncoly+1)/2)){
        realized[, k] <- apply.quarterly(cross1[,k], sum)
      }
      chol2 <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nquarter)
      for (j in 1:nmonth){
        rcovmat <- invvech(as.matrix(realized[j,]))
        chol1 <- t(chol(rcovmat))
        chol2[j,] <- vech(chol1)
      }
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
        realized[, k] <- apply.quarterly(cross1[,k], sum)
      }
      chol2 <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nyear)
      for (j in 1:nyear){
        rcovmat <- invvech(as.matrix(realized[j,]))
        chol1 <- t(chol(rcovmat))
        chol2[j,] <- vech(chol1)
      }
    }}

names1 <- to('y', (ncoly*(ncoly+1)/2), same.size = F)

colnames(realized) <- names1
colnames(chol2) <- names1

if(cholesky == T){
  results <- list(realized, chol2)
  names(results) <- c('Realized Covariances', 'Cholesky Factors')
}else{
  results <- list(realized)
  names(results) <- 'Realized Covariances'
}
return(results)
}
