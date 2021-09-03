#' Log-Likelihood method
#'
#' Returns the log-Likelihood of a VLSTAR object.
#'
#'   The log-likelihood of a VLSTAR model is defined as: \deqn{\log l(y_t|I_t;\theta)=-\frac{T\tilde{n}}{2}\ln(2\pi)-\frac{T}{2}\ln|\Omega|-\frac{1}{2}\sum_{t=1}^{T}(y_t-\tilde{G}_tB\,z_t)'\Omega^{-1}(y_t-\tilde{G}_tB\,z_t)}
#'
#'@usage \method{logLik}{VLSTAR}(object, type = c('Univariate', 'Multivariate'), \dots)
#'@author Andrea Bucci
#'@param object An object of class \sQuote{\code{VLSTAR}} obtained through \command{VLSTAR()}.
#'@param type Type of Log-Likelihood to be showed (univariate or multivariate).
#'@param \dots further arguments to be passed to and from other methods
#'@return An object with class attribute \code{logLik}.
#'@seealso \code{\link{VLSTAR}}
#'@aliases logLik logLik.VLSTAR
#'@export


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
