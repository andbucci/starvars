dummyfilter <- function(y, date1, date2){
  if(is.null(date1) & is.null(date2)){
    stop('Please provide a vector of single or consecutive dates')}
  if(dim(y)[2] > 1){
    stop('Function allows only univariate time series')
  }
  if(is(y, 'xts') == F){
    stop('Data should be of class "xts".')
  }
  data <- as.matrix(y)
  nrowy <- nrow(data)
  ndate1 <- length(date1) #single dates
  ndate2 <- length(date2) #consecutive dates
  dummy1 <- matrix(0, nrow= nrowy, ncol = (ndate1+ndate2))
  dummy1 <- zoo(dummy1, index(y))
  if(!is.null(date1) & !is.null(date2)){
  for(i in 1:ndate1){
   dummy1[index(dummy1)==date1[i], i] <- 1 
  }
  for(j in 1:ndate2){
    dummy1[index(dummy1)==date1[i], (ndate2+j)] <- 1
  } else if(!is.null(date1) & is.null(date2)){
    for(i in 1:ndate1){
      dummy1[index(dummy1)==date1[i], i] <- 1 
    }
  } else if(is.null(date1) & !is.null(date2)){
  lm1 <- lm(dato ~ dummy1)
  res1 <- residuals(lm1)
  res1
}