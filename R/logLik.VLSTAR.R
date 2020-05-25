logLik.VLSTAR <- function(object, type = c('Univariate', 'Multivariate'), ...){
  if(!type %in% c('Univariate', 'Multivariate'))
    stop('Please, provide a valid argument')
  df = (nrow(object$yoriginal)*object$m-object$m*(ncol(object$Data[[2]])))
  obs <- nrow(object$Data[[2]])
  if(type == 'Univariate'){
    r <- object$LL
  }else{
    r <- object$MultiLL
  }
  class(r) <- "logLik"
  attr(r, "df") <- df
  attr(r, "nobs") <- obs
  return(r)
}
