#' @title Estimation of Coefficients in a Vector Smooth Transition Autoregressive Model with m regimes through maximum likelihood
#'
#' @description This package  allows to estimate the coefficient in a VLSTAR model.
#'
#'
#' @examples  VLSTAR.lm(data, p = 1, st = yy)
#'

VLSTAR.lm <- function(y1, x1 = NULL, p = NULL,
                      m = 2, st = NULL, constant = T,
                      n.combi = 50, n.iter = 500,
                      starting = NULL, epsilon = 10^(-3),
                      exo = T){
  y <- as.matrix(y1)
  if (any(is.na(y)))
    stop("\nNAs in y.\n")
  if(m < 2)
    stop('The number of regimes should be greater than one.')
  if(is.null(st))
    stop('The transition variable must be supplied.')
  if(length(y[,1]) != length(as.matrix(x1[,1])) | length(st) != length(as.matrix(x1[,1])) | length(y[,1]) != length(st))
    stop('The length of the variables does not match!')
  if(is.null(p) | p < 1){
    stop('Please, specify a valid lag order.')
  }
  if (ncol(y) < 2)
    stop("The matrix 'y' should contain at least two variables. For univariate analysis consider lstar() function in this package.\n")
  if (is.null(colnames(y))) {
    colnames(y) <- paste("y", 1:ncol(y), sep = "")
    warning(paste("No column names supplied in y, using:",
                  paste(colnames(y), collapse = ", "), ", instead.\n"))
  }
  colnames(y) <- make.names(colnames(y))
  ##Definition of dimensions, creating variable x with constant
  yt <- zoo(y)
  ylag <- stats::lag(yt, -c(1:p))
  ylag <- as.matrix(ylag)
  y <- y[-p, ]
  ncoly <- ncol(y1)
  ncolylag <- ncol(ylag)
  nrowy <- nrow(y1)
  ncolx1 <- ncol(x1)
  const <- rep(1, (nrowy-p))
  if (constant == T){
    if(exo == T){
      x <- as.matrix(cbind(ylag, constant, x1[-p,]))
    }else{
      x <- as.matrix(cbind(ylag, constant))
    }
    ncolx <- ncol(x)
  }  else{
    x <- as.matrix(x1)
  }
  nrowx <- nrow(x)
  ncolx <- ncol(x)
  param.init <- list()
  param.init$gamma <- rep(1L, ncoly)
  param.init$cg <- colMeans(y)
  q <- ncol(x)-ncolylag

  if (is.null(starting)){
    COMBI <- list()
    GAMMA <- list()
    CJ <- list()
  for(t in 1:(m-1)){
  #Starting values for c and gamma
  sumsq <- matrix(ncol = ncoly, nrow = (n.combi*n.combi))
  gamma <- matrix(ncol = ncoly, nrow = n.combi)
  cj <- matrix(ncol = ncoly, nrow = n.combi)
  gamma[1, ] <- param.init$gamma
  cj[1,] <- param.init$c
  GAMMA[[t]] <- gamma*(1+rnorm(1, mean = 0, sd = 5))
  CJ[[t]] <- cj*(1+rnorm(1, mean = 0, sd = 5))
  #Grid for c and gamma
  rangey <- matrix(nrow = ncoly, ncol = 2)
  combi <- list()
  for(j in 1:ncoly){
    rangey[j,] <- range(y[,j])
    CJ[[t]][2:n.combi,j] <- seq(from = rangey[j,1]*1.1, to = rangey[j,2], length.out = (n.combi-1))
    GAMMA[[t]][2:n.combi, j] <- seq(from = 0L, to = 100L, length.out = (n.combi-1))
    combi[[j]] <- expand.grid(CJ[[t]][,j], GAMMA[[t]][,j])
  }
  COMBI[[t]] <- combi
  }

  #NLS for each combination of c and gamma
  ssq <- matrix(ncol =  ncoly, nrow = n.combi*n.combi)
  coeff <- matrix(ncol = (ncolx*m*ncoly), nrow = n.combi*n.combi)

  for (l in 1:(n.combi*n.combi)){
    tryCatch({
      In <- diag(ncoly)
      glog <- matrix(ncol=ncoly, nrow = nrowy)
      GT <- list()
      Gtilde <- list()
      GG <- list()
      XX <- list()
      GGXX <- list()
      XY <- list()
      XYG <- list()
      kro <- list()
      ggxx <- matrix(ncol = (q+ncoly)*m*ncoly, nrow = (q+ncoly)*m*ncoly)
      for (i in 1:nrowx){
        for(t in 1:(m-1)){
        for (j in 1 : ncoly){
          glog[i,j] <- (1+exp(-COMBI[[t]][[j]][l,2]*(st[i]-COMBI[[t]][[j]][l,1])))^(-1)
        }
        Gt <- diag(glog[i,])
        GT[[t]] <- Gt
      }
      Gtilde[[i]] <- t(cbind(In, list.cbind(GT)))
        GG[[i]] <- Gtilde[[i]]%*%t(Gtilde[[i]])
        XX[[i]] <- x[i,] %*%t(x[i,])
        GGXX[[i]] <- kronecker(GG[[i]], XX[[i]])
        XY[[i]] <- x[i, ]%*%t(y[i,])
        XYG[[i]] <- vec(XY[[i]]%*%t(Gtilde[[i]]))
        kro[[i]] <- kronecker(Gtilde[[i]], x[i,])
      }
      ggxx <- Reduce(`+`, GGXX)/nrowx
      M <- t(do.call("cbind", kro))
      Y <- vec(t(y))
      Bhat <- ginv(t(M)%*%M)%*%t(M)%*%Y
      BB <- invvec(Bhat, ncol = (m*ncoly), nrow = (ncoly + q))
      resi <- list()
      Ehat <- matrix(NA, ncol = ncoly, nrow = nrowy)
      for (o in 1:nrowx){
        resi[[o]] <- y[o, ] - t(Gtilde[[o]])%*%t(BB)%*%x[o,]
      }
      Ehat <- t(do.call("cbind", resi))
      SSQ <- colSums(Ehat^2)
      ssq[l,] <- SSQ
      coeff[l, ] <- Bhat
      print(l)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }

  #c and gamma minimazing the sum of squared residuals for each equation

  cgamma <- matrix(nrow = ncoly, ncol = 2)
  CGAMMA <- list()
  for (t in 1:(m-1)){
    for(j in 1:ncoly){
      cgamma[j,] <- as.matrix(COMBI[[t]][[j]][which.min(ssq[,j]),])
    }
    CGAMMA[[t]] <- cgamma
  }


  #Definition of c0 gamma0 for maximum likelihood
  PARAM <- list()
  for (t in 1:(m-1)){
    cj <- as.matrix(CGAMMA[[t]][,1])
    gamma <- as.matrix(CGAMMA[[t]][,2])
    PARAM[[t]] <- cbind(gamma, cj)
  }

}else{
  PARAM <- starting
}

  #NLS Estimation of Bhat and Omegahat to be used in the first iteration of maximum likelihood
  glog <- matrix(ncol=ncoly, nrow = nrowy)
  GT <- list()
  Gtilde <- list()
  GG <- list()
  XX <- list()
  GGXX <- list()
  XY <- list()
  XYG <- list()
  kro <- list()
  ggxx <- matrix(ncol = (q+ncoly)*m*ncoly, nrow = (q+ncoly)*m*ncoly)
  for (i in 1:nrowx){
    for(t in 1:(m-1)){
    for (j in 1 : ncoly){
      (1L+exp(-PARAM[[t]][j,1]*(st[i]-PARAM[[t]][j,2])))^(-1)
    }
    Gt <- diag(glog[i,])
    GT[[t]] <- Gt
  }
  Gtilde[[i]] <- t(cbind(In, list.cbind(GT)))
    GG[[i]] <- Gtilde[[i]]%*%t(Gtilde[[i]])
    XX[[i]] <- x[i,] %*%t(x[i,])
    GGXX[[i]] <- kronecker(GG[[i]], XX[[i]])
    XY[[i]] <- x[i, ]%*%t(y[i,])
    XYG[[i]] <- vec(XY[[i]]%*%t(Gtilde[[i]]))
    kro[[i]] <- kronecker(Gtilde[[i]], x[i,])
  }
  ggxx <- Reduce(`+`, GGXX)/nrowx
  M <- t(do.call("cbind", kro))
  Y <- vec(t(y))
  #Estimated coefficients
  Bhat <- ginv(t(M)%*%M)%*%t(M)%*%Y
  BB <- invvec(Bhat, ncol = (m*ncoly), nrow = (ncoly + q))
  resi <- list()
  Ehat1 <- matrix(NA, ncol = ncoly, nrow = nrowy)
  for (i in 1:nrowx){
    resi[[i]] <- y[i, ] - t(Gtilde[[i]])%*%t(BB)%*%x[i,]
  }
  Ehat1 <- t(do.call("cbind", resi))
  #Estimated covariance matrix
  Omegahat <- (t(Ehat1)%*%Ehat1)/(nrowy)

  glog <- matrix(ncol=ncoly, nrow = nrowy)
  GT <- list()
  Gtilde <- list()
  GG <- list()
  XX <- list()
  XY <- list()
  XYOG <- list()
  kro <- list()
  PsiOmegaPsi <- list()
  ggxx <- matrix(ncol = (q+ncoly)*m*ncoly, nrow = (q+ncoly)*m*ncoly)
  for (i in 1:nrowx){
    for(t in 1:(m-1)){
    for (j in 1 : ncoly){
      glog[i,j] <- (1+exp(-cgam[j,1]*(st[i]-cgam[j,2])))^(-1)
    }
    Gt <- diag(glog[i,])
    GT[[t]] <- Gt
  }
  Gtilde[[i]] <- t(cbind(In, list.cbind(GT)))
    GG[[i]] <- Gtilde[[i]]%*%t(Gtilde[[i]])
    XX[[i]] <- x[i,] %*%t(x[i,])
    XY[[i]] <- x[i, ]%*%t(y[i,])
    XYOG[[i]] <- vec(XY[[i]]%*%ginv(Omegahat)%*%t(Gtilde[[i]]))
    PsiOmegaPsi[[i]] <- Gtilde[[i]]%*%ginv(Omegahat)%*%t(Gtilde[[i]])
    kro[[i]] <- kronecker(PsiOmegaPsi[[i]], XX[[i]])
  }
  xyog <- Reduce(`+`, XYOG)/nrowy
  kroxx <- Reduce(`+`, kro)/nrowy
  Bhat <- t(xyog%*%t(kroxx))
  BB <- invvec(Bhat, ncol = (m*ncoly), nrow = (ncoly + q))

  #BB0 <- BB
  #Omega0 <- Omegahat

  #Log-likelihood to be optimized
  loglike <- function(param){
    gamma <- param[1:((m-1)*ncoly)]
    c <- param[(ncoly*(m-1)+1):length(param)]
    gamma1 <- matrix(gamma, ncol = (m-1))
    c1 <- matrix(c, ncol = (m-1))
    glog <- matrix(ncol=ncoly, nrow = nrowy)
    GT <- list()
    Gtilde <- list()
    dify <- matrix(ncol = 1, nrow = nrow(y))
    ncoly <- dim(y)[2]
    In <- diag(ncoly)
    for (z in 1:nrow(y)){
      for(t in 1:(m-1)){
        for (o in 1:ncoly){
          gammao <- gamma[o]
          co <- c[o]
          glog[z,o] <- (1L+exp(-gammao*(st[z]-co)))^(-1)}
        Gt <- diag(glog[z,])
        GT[[t]] <- Gt
      }
      Gtilde[[z]] <- t(cbind(In, list.cbind(GT)))
      dify[z] <-  t(y[z, ] - t(Gtilde[[z]])%*%t(BB)%*%x[z,])%*%ginv(Omegahat)%*%(y[z, ] - t(Gtilde[[z]])%*%t(BB)%*%x[z,])
    }
    sumdif <- sum(dify)
    logll <- -(nrowy*log(det(Omegahat))/2L) - sumdif/2L  - (nrowy*ncoly/2L)*log(2L*pi)
    return(-logll)
  }

  #Maximum likelihood calculation
  loglike1 <- function(BB, param, Omega){
    gamma <- param[1:((m-1)*ncoly)]
    c <- param[(ncoly*(m-1)+1):length(param)]
    gamma1 <- matrix(gamma, ncol = (m-1))
    c1 <- matrix(c, ncol = (m-1))
    glog <- matrix(ncol=ncoly, nrow = nrowy)
    GT <- list()
    Gtilde <- list()
    dify <- matrix(ncol = 1, nrow = nrow(y))
    ncoly <- dim(y)[2]
    In <- diag(ncoly)
    for (z in 1:nrow(y)){
      for(t in 1:(m-1)){
        for (o in 1:ncoly){
          gammao <- gamma[o]
          co <- c[o]
          glog[z,o] <- (1L+exp(-gammao*(st[z]-co)))^(-1)}
        Gt <- diag(glog[z,])
        GT[[t]] <- Gt
      }
      Gtilde[[z]] <- t(cbind(In, list.cbind(GT)))
      dify[z] <-  t(y[z, ] - t(Gtilde[[z]])%*%t(BB)%*%x[z,])%*%ginv(Omegahat)%*%(y[z, ] - t(Gtilde[[z]])%*%t(BB)%*%x[z,])
    }
    sumdif <- sum(dify)
    logll <- -(nrowy*log(det(Omegahat))/2L) - sumdif/2L  - (nrowy*ncoly/2L)*log(2L*pi)
    return(-logll)
  }

  #Inizialization of iter
  iter <- 0
  ll0 <- 10^(4)
  epsi <- epsilon
  err <- list()
  loglik1 <- NULL
  bbhat <- list()
  omega <- list()
  errdif <- 10^5
  PARAM1 <- list()
  for(t in 1:(m-1)){
    PARAM1[[t]] <- as.data.frame(PARAM[[t]])
  }
  param <- rbindlist(PARAM1)
  param <- vec(param)

  cat('Maximum likelihood estimation\n')

  #Convergence algorithm
  #1. Maximum likelihood estimation of gamma and c with NLS estimates of Bhat and Omegahat
  #2. Maximum likelihood estimation of Bhat with new values of gamma and c
  #3. Convergence check
  #4. 1-2-3 until convergence
  while (iter < n.iter & errdif > epsi){
    Sys.sleep(0)
    iter <- iter+1
    #Parameters
    low1 <- replicate(ncoly, 0)
    #1.Maximum likelihood estimation of gamma and c
    param1 <- optim(par = param, fn = loglike, lower = c(low1, apply(y, 2, min)),
                    method="L-BFGS-B")

    cgam1 <- matrix(param1$par, ncol = 2, nrow = ncoly)

    #2.Maximum likelihood estimation of Bhat with new values of gamma and c
    glog <- matrix(ncol=ncoly, nrow = nrowy)
    GT <- list()
    Gtilde <- list()
    GG <- list()
    XX <- list()
    XY <- list()
    XYOG <- list()
    kro <- list()
    PsiOmegaPsi <- list()
    ggxx <- matrix(ncol = (q+ncoly)*m*ncoly, nrow = (q+ncoly)*m*ncoly)
    for (i in 1:nrowx){
      for(t in 1:(m-1)){
      for (j in 1 : ncoly){
        glog[i,j] <- (1+exp(-cgam1[j,1]*(st[i]-cgam1[j,2])))^(-1)
      }
        Gt <- diag(glog[i,])
        GT[[t]] <- Gt
      }
      Gtilde[[i]] <- t(cbind(In, list.cbind(GT)))
      GG[[i]] <- Gtilde[[i]]%*%t(Gtilde[[i]])
      XX[[i]] <- x[i,] %*%t(x[i,])
      XY[[i]] <- x[i, ]%*%t(y[i,])
      XYOG[[i]] <- vec(XY[[i]]%*%ginv(Omegahat)%*%t(Gtilde[[i]]))
      PsiOmegaPsi[[i]] <- Gtilde[[i]]%*%ginv(Omegahat)%*%t(Gtilde[[i]])
      kro[[i]] <- kronecker(PsiOmegaPsi[[i]], XX[[i]])
    }
    xyog <- Reduce(`+`, XYOG)/nrowy
    kroxx <- Reduce(`+`, kro)/nrowy
    Bhat <- t(xyog%*%t(kroxx))
    BB <- invvec(Bhat, ncol = (m*ncoly), nrow = (ncoly + q))
    resi <- list()
    Ehat1 <- matrix(NA, ncol = ncoly, nrow = nrowy)
    for (i in 1:nrowx){
      resi[[i]] <- y[i, ] - t(Gtilde[[i]])%*%t(BB)%*%x[i,]
    }
    Ehat1 <- t(do.call("cbind", resi))
    Omegahat <- (t(Ehat1)%*%Ehat1)/(nrowy)

    #3. Convergence check
    ll1 <- loglike1(BB, param1$par, Omega = Omegahat)
    if((ll1-ll0)>0){
      param <- as.matrix(vec(cgam1))
    }
    err <- abs(ll1 - ll0)
    ll0 <- ll1
    if(iter ==1) {
      plot(ll0~1, xlim = c(0, n.iter), ylim = c(0, (ll0+ll0*0.1)))
    }
    points(ll0~iter)
    loglik1[iter] <- ll1
    bbhat[[iter]] <- BB
    omega[[iter]] <- Omegahat

    cat(paste("iteration", iter, "complete\n"))

    print(ll0)

    if (err<epsi | iter == n.iter) cat('Converged\n')}

  #Bhat1 <- bbhat[[which.min(loglik1)]]
  residuals1 <- t(do.call("cbind", resi))
  varhat <- diag(omega[[iter]])
  bb1 <- bbhat[[iter]][,1:ncoly]
  bb2 <- list()
  for(t in 1:(m-1)){
    bb2[[t]] <- as.data.frame(bbhat[[iter]][,(ncoly*(t-1)+1):(ncoly*t)] +
                                bbhat[[iter]][,(ncoly*(t)+1):(ncoly*(t+1))])
  }
  bb4 <- as.matrix(rbindlist(bb2))
  BBhat <- rbind(bb1, bb4)
  colnames(cgam1) <- c('gamma', 'c')
  covbb <- matrix(nrow = m*ncolx, ncol = ncoly)
  ttest <- matrix(nrow = m*ncolx, ncol = ncoly)
  pval <- matrix(nrow = m*ncolx, ncol = ncoly)
  ee <- matrix(nrow = m*ncolx, ncol = ncoly)
  for (j in 1 : ncoly){
    #covbb[,j] <- diag(ginv(t(x[[j]])%*%XX[[j]])*sqrt(varhat[j]))
    covbb[,j] <- sqrt(diag(ginv(t(x) %*%x))*varhat[j])
    ttest[,j] <- BBhat[,j]/covbb[,j]
    pval[,j] <- 2*pt(abs(ttest[,j]),df=(nrowy*m-m*(ncolx)), lower = F)
  }
  signifi <- matrix('',nrow = nrow(pval), ncol = ncol(pval))
  for(i in 1:nrow(pval)){
    for(j in 1:ncol(pval)){
      if(pval[i,j] < 0.01)
      {signifi[i,j] <- '***'}
      else if(pval[i,j] <=0.05 & pval[i,j]>0.01){
        signifi[i,j] <- '**'
      }
      else if(pval[i,j] <= 0.1 & pval[i,j]>0.05){
        signifi[i,j] <- '*'
      }
    }
  }

  bhat1 <- matrix(nrow = nrow(pval), ncol = ncol(pval))
  for(i in 1:nrow(pval)){
    for(j in 1:ncol(pval)){
      bhat1[i,j] <- paste(round(BBhat[i,j],4), signifi[i,j], "&", "(", round(covbb[i,j],4),
                          ")", sep = '')
    }
  }
  loglike2 <- function(y, resid1, omega){
    nrowy <- nrow(as.matrix(y))
    logll <- -(nrowy/2)*log(2*pi) -(nrowy/2)*log(omega) - (t(resid1)%*%resid1)/(2*omega)
    return(logll)
  }
  residui <- residuals1
  omega1 <- diag(omega[[iter]])
  k <- nrow(BBhat)
  ll2 <- NULL
  AIC1 <- NULL
  BIC1 <- NULL
  for (l in 1:ncoly){
    ll2[l] <- loglike2(y1[,l], residui[,l], omega1[l])
    AIC1[l] <- 2*k - 2*ll2[l]
    BIC1[l] <- -2*ll2[l] + k*log(nrowy)
  }
  names <- list()
  for(j in 1:m){
    names1[[j]] <- as.data.frame(paste(colnames(x), 'm_', j))
  }
  names1 <- as.matrix(rbindlist(tmp))
  rownames(BBhat) <- names1
  colnames(BBhat) <- colnames(y)
  rownames(bhat1) <- names1
  colnames(bhat1) <- colnames(y)
  modeldata <- list(y, x)
  results <- list(BBhat, covbb, ttest, pval, cgam1, omega[[iter]], fitte, residuals1, bhat1, ll1, ll2, AIC1, BIC1, Gt, modeldata, BB)
  names(results) <- c('Bhat','St.Dev.', 't-test', 'pval', 'C-gamma', 'Omega', 'Fitted', 'Residuals', 'Output', 'Multivariate Log-Likelihood', 'Log-Likelihood', 'AIC', 'BIC', 'Gtilde', 'Data', 'B')
  return(results)
  }
