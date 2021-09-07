#'Multivariate log-likelihood
#'
#'Log-likelihood to be optimized in the ML method and used to check convergence in both methods
#'@keywords internal
#'
loglike <- function(param, data){
  y = data$y
  Omegahat = data$Omegahat
  gamma <- param[1:((data$m-1)*ncol(y))]
  c <- param[(ncol(y)*(data$m-1)+1):length(param)]
  glog <- rep(0, ncol(y))
  GT <- list()
  dify <- matrix(ncol = 1, nrow = nrow(y))
  for (z in 1:nrow(y)){
    for(t in 1:(data$m-1)){
      for (o in 1:ncol(y)){
        glog[o] <- (1L+exp(-gamma1[o,t]*(data$st[z]-c1[o,t])))^(-1)}
      ifelse(data$singlecgamma == TRUE, GT[[t]] <- diag(rep(glog[1], ncol(y))), GT[[t]] <- diag(glog))
    }
    Gtilde <- t(cbind(diag(ncol(y)), do.call(cbind,GT)))
    dify[z] <-  t(y[z, ] - t(Gtilde)%*%t(data$BB)%*%data$x[z,])%*%MASS::ginv(Omegahat)%*%(y[z, ] - t(Gtilde)%*%t(data$BB)%*%data$x[z,])
  }
  sumdif <- sum(dify)
  logll <- -(nrow(y)*log(det(Omegahat))/2L) - sumdif/2L  - (nrow(y)*ncol(y)/2L)*log(2L*pi)##Normal distribution assumed
  return(-logll)
}

#'Sum of squared error
#'
#'Sum of squared error to be optimized in the NLS method
#'@keywords internal
#'
SSQ <- function(param, data){
  y <- data$y
  gamma <- param[1:((data$m-1)*ncol(y))]
  c <- param[(ncol(y)*(data$m-1)+1):length(param)]
  gamma1 <- matrix(gamma, ncol = (data$m-1))
  c1 <- matrix(c, ncol = (data$m-1))
  glog <- rep(0, ncol(y))
  GT <- list()
  dify <- matrix(ncol = 1, nrow = nrow(y))
  for (z in 1:nrow(y)){
    for(t in 1:(data$m-1)){
      for (o in 1:ncol(y)){
        glog[o] <- (1L+exp(-gamma1[o,t]*(data$st[z]-c1[o,t])))^(-1)}
      ifelse(data$singlecgamma == TRUE, GT[[t]] <- diag(rep(glog[1], ncol(y))), GT[[t]] <- diag(glog))
    }
    #Gtilde <- t(cbind(diag(ncol(y)), do.call(cbind,GT)))
    dify[z] <-  t(y[z, ] - t(t(cbind(diag(ncol(y)), do.call(cbind,GT))))%*%t(data$BB)%*%data$x[z,])%*%(y[z, ] - t(t(cbind(diag(ncol(y)), do.call(cbind,GT))))%*%t(data$BB)%*%data$x[z,])
  }
  sumdif <- sum(dify)
  return(sumdif)
}
