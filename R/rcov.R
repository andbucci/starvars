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
  year(end_date)-
    year(start_date)+1
}
nmonth <- elapsed_months(end_date, start_date)
nquarter <- elapsed_quarters(end_date, start_date)
nyear <- elapsed_years(end_date, start_date)


retu <- matrix(ncol = ncoly, nrow = nrowy)
if(make.ret == TRUE){
    if(freq == 'daily'){
      realized <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nday)
      for (k in 1:nday){
        realized[k,] <- vech(rCov(rdata = data[as.character(days[k])], makeReturns = TRUE))
      }
      chol2 <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nday)
      for (j in 1:nday){
        rcovmat <- invvech(as.matrix(realized[j,]))
        chol1 <- t(chol(rcovmat))
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
        chol1 <- t(chol(rcovmat, pivot = T))
        chol2[j,] <- vech(chol1)
      }
      realized <- as.xts(realized, order.by = seq(start_date, end_date, by = 'quarter'))
      chol2 <- as.xts(chol2, order.by = seq(start_date, end_date, by = 'quarter'))
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
      realized <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nday)
      for (k in 1:nday){
        realized[k,] <- vech(rCov(rdata = data[as.character(days[k])], makeReturns = FALSE))
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
        chol1 <- t(chol(rcovmat, pivot = T))
        chol2[j,] <- vech(chol1)
      }
      realized <- as.xts(realized, order.by = seq(start_date, end_date, by = 'month'))
      chol2 <- as.xts(chol2, order.by = seq(start_date, end_date, by = 'month'))
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

if(cholesky == TRUE){
  results <- list(realized, chol2)
  names(results) <- c('Realized Covariances', 'Cholesky Factors')
}else{
  results <- list(realized)
  names(results) <- 'Realized Covariances'
}
return(results)
}
