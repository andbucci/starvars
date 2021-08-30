logLik.VLSTAR <- function(object, type = c('Univariate', 'Multivariate'), ...){
  type <- match.arg(type)
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
