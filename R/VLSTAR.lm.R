#' @title Estimation of Coefficients in a Vector Smooth Transition Autoregressive Model with 2 regimes through maximum likelihood
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
  require(vars)
  require(nloptr)
  require(matrixcalc)
  require(MASS)
  library(ks)
  require(maxLik)
  library(stats4)
  y <- as.matrix(y1)
  if (any(is.na(y)))
    stop("\nNAs in y.\n")
  if(m < 2)
    stop('The number of regimes should be greater than one.')
  if(is.null(st))
    stop('The transition variable must be supplied.')
  if(length(y[,1]) != length(as.matrix(x1[,1])) | length(st) != length(as.matrix(x1[,1])) | length(y[,1]) != length(st))
    stop('The length of the variables does not match!')
  if(is.null(p)){
    warning('A single lag will be used!')
    p <- 1
  }
  if (ncol(y) < 2)
    stop("The matrix 'y' should contain at least two variables. For univariate analysis consider lstar() function in this package.\n")
  if (is.null(colnames(y))) {
    colnames(y) <- paste("y", 1:ncol(y), sep = "")
    warning(paste("No column names supplied in y, using:",
                  paste(colnames(y), collapse = ", "), ", instead.\n"))
  }

  ##Definition of dimensions, creating variable x with constant
  ncoly <- dim(y1)[2]
  nrowy <- dim(y1)[1]
  ncolx1 <- dim(x1)[2]
  const <- rep(1, nrowy)
  if (constant == T){
    x <- as.matrix(cbind(x1[,1:ncoly], const, x1[,(ncoly+1):ncolx1]))
    ncolx <- dim(x)[2]
    #xvar <- x[,(ncoly+2):ncolx]
  }  else{
    x <- as.matrix(x1)
    xvar <- x[,(ncoly+1):ncolx]}
  nrowx <- dim(x)[1]
  ncolx <- dim(x)[2]
  param.init <- list()
  param.init$gamma <- rep(1, ncoly)
  param.init$cg <- colMeans(y)
  q <- dim(x)[2]-ncoly

  if (is.null(starting)){
  #Starting values for c and gamma
  sumsq <- matrix(ncol = ncoly, nrow = (n.combi*n.combi))
  gamma <- matrix(ncol = ncoly, nrow = n.combi)
  cj <- matrix(ncol = ncoly, nrow = n.combi)
  gamma[1, ] <- param.init$gamma
  cj[1,] <- param.init$c
  #Grid for c and gamma
  rangey <- matrix(nrow = ncoly, ncol = 2)
  combi <- list()
  for(j in 1:ncoly){
    rangey[j,] <- range(y[,j])
    cj[2:n.combi,j] <- seq(from = rangey[j,1]*1.1, to = rangey[j,2], length.out = (n.combi-1))
    gamma[2:n.combi, j] <- seq(from = 0, to = 100, length.out = (n.combi-1))
    combi[[j]] <- expand.grid(cj[,j], gamma[,j])
  }

  #NLS for each combination of c and gamma
  ssq <- matrix(ncol =  ncoly, nrow = n.combi*n.combi)
  coeff <- matrix(ncol = (ncolx*m*ncoly), nrow = n.combi*n.combi)

  for (l in 1:(n.combi*n.combi)){
    tryCatch({
      In <- diag(ncoly)
      glog <- matrix(ncol=ncoly, nrow = nrowy)
      Gt <- list()
      Gtilde <- list()
      GG <- list()
      XX <- list()
      GGXX <- list()
      XY <- list()
      XYG <- list()
      kro <- list()
      ggxx <- matrix(ncol = (q+ncoly)*m*ncoly, nrow = (q+ncoly)*m*ncoly)
      for (i in 1:nrowx){
        for (j in 1 : ncoly){
          glog[i,j] <- (1+exp(-combi[[j]][l,2]*(st[i]-combi[[j]][l,1])))^(-1)
        }
        Gt[[i]] <- diag(glog[i,])
        Gtilde[[i]] <- t(cbind(In, Gt[[i]]))
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
  for(j in 1:ncoly){
    cgamma[j,] <- as.matrix(combi[[j]][which.min(ssq[,j]),])
  }

  cgam <- cgamma

  #Definition of c0 gamma0 for maximum likelihood
  cj <- as.matrix(cgam[,1])
  gamma <- as.matrix(cgam[,2])
  param <- rbind(gamma, cj)}else{
    param <- starting
  }

  #NLS Estimation of Bhat and Omegahat to be used in the first iteration of maximum likelihood
  glog <- matrix(ncol=ncoly, nrow = nrowy)
  Gt <- list()
  Gtilde <- list()
  GG <- list()
  XX <- list()
  GGXX <- list()
  XY <- list()
  XYG <- list()
  kro <- list()
  ggxx <- matrix(ncol = (q+ncoly)*m*ncoly, nrow = (q+ncoly)*m*ncoly)
  for (i in 1:nrowx){
    for (j in 1 : ncoly){
      glog[i,j] <- (1L+exp(-cgam[j,2]*(st[i]-cgam[j,1])))^(-1)
    }
    Gt[[i]] <- diag(glog[i,])
    Gtilde[[i]] <- t(cbind(In, Gt[[i]]))
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
  Gt <- list()
  Gtilde <- list()
  GG <- list()
  XX <- list()
  XY <- list()
  XYOG <- list()
  kro <- list()
  PsiOmegaPsi <- list()
  ggxx <- matrix(ncol = (q+ncoly)*m*ncoly, nrow = (q+ncoly)*m*ncoly)
  for (i in 1:nrowx){
    for (j in 1 : ncoly){
      glog[i,j] <- (1+exp(-cgam[j,1]*(st[i]-cgam[j,2])))^(-1)
    }
    Gt[[i]] <- diag(glog[i,])
    Gtilde[[i]] <- t(cbind(In, Gt[[i]]))
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
    gamma <- param[1:ncoly]
    c <- param[(ncoly+1):(ncoly*2)]
    glog <- matrix(ncol=ncoly, nrow = nrowy)
    Gt <- list()
    Gtilde <- list()
    dify <- matrix(ncol = 1, nrow = nrow(y))
    ncoly <- dim(y)[2]
    In <- diag(ncoly)
    for (z in 1:nrow(y)){
      for (o in 1:ncoly){
        gammao <- gamma[o]
        co <- c[o]
        glog[z,o] <- (1L+exp(-gammao*(st[z]-co)))^(-1)}
      Gt[[z]] <- diag(glog[z,])
      Gtilde[[z]] <- t(cbind(In, Gt[[z]]))
      dify[z] <-  t(y[z, ] - t(Gtilde[[z]])%*%t(BB)%*%x[z,])%*%ginv(Omegahat)%*%(y[z, ] - t(Gtilde[[z]])%*%t(BB)%*%x[z,])
    }
    sumdif <- sum(dify)
    logll <- -(nrowy*log(det(Omegahat))/2) - sumdif/2  - (nrowy*ncoly/2)*log(2*pi)
    return(-logll)
  }

  #Maximum likelihood calculation
  loglike1 <- function(BB, param, Omega){
    gamma <- param[1:ncoly]
    c <- param[(ncoly+1):(ncoly*2)]
    glog <- matrix(ncol=ncoly, nrow = nrowy)
    Gt <- list()
    Gtilde <- list()
    dify <- matrix(ncol = 1, nrow = nrow(y))
    ncoly <- dim(y)[2]
    In <- diag(ncoly)
    for (z in 1:nrow(y)){
      for (o in 1:ncoly){
        gammao <- gamma[o]
        co <- c[o]
        glog[z,o] <- (1L+exp(-gammao*(st[z]-co)))^(-1)}
      Gt[[z]] <- diag(glog[z,])
      Gtilde[[z]] <- t(cbind(In, Gt[[z]]))
      dify[z] <-  t(y[z, ] - t(Gtilde[[z]])%*%t(BB)%*%x[z,])%*%ginv(Omegahat)%*%(y[z, ] - t(Gtilde[[z]])%*%t(BB)%*%x[z,])
    }
    sumdif <- sum(dify)
    logll <- -(nrowy*log(det(Omegahat))/2) - sumdif/2  - (nrowy*ncoly/2)*log(2*pi)
    return(-logll)
  }

  #Inizialization of iter
  iter <- 0
  ll0 <- 10^(4)
  epsi <- 10^(-3)
  err <- list()
  loglik1 <- NULL
  bbhat <- list()
  omega <- list()
  errdif <- 10^5
  param <- vec(cgam)

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
    Gt <- list()
    Gtilde <- list()
    GG <- list()
    XX <- list()
    XY <- list()
    XYOG <- list()
    kro <- list()
    PsiOmegaPsi <- list()
    ggxx <- matrix(ncol = (q+ncoly)*m*ncoly, nrow = (q+ncoly)*m*ncoly)
    for (i in 1:nrowx){
      for (j in 1 : ncoly){
        glog[i,j] <- (1+exp(-cgam1[j,1]*(st[i]-cgam1[j,2])))^(-1)
      }
      Gt[[i]] <- diag(glog[i,])
      Gtilde[[i]] <- t(cbind(In, Gt[[i]]))
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
  #Minimum of minus log likelihood
  Bhat1 <- bbhat[[which.min(loglik1)]]
  residuals1 <- t(do.call("cbind", resi))
  varhat <- diag(omega[[which.min(loglik1)]])
  Bhat2 <- c(Bhat1[1:(ncolx*m*ncoly/2)], (Bhat1[1:(ncolx*m*ncoly/2)]+Bhat1[((ncolx*m*ncoly/2)+1):(ncolx*m*ncoly)]))
  BBhat <- invvec(Bhat2, ncol = ncoly, nrow = (ncolx*m))
  covbb <- matrix(nrow = m*ncolx, ncol = ncoly)
  ttest <- matrix(nrow = m*ncolx, ncol = ncoly)
  pval <- matrix(nrow = m*ncolx, ncol = ncoly)
  ee <- matrix(nrow = m*ncolx, ncol = ncoly)
  for (j in 1 : ncoly){
    covbb[,j] <- sqrt(diag(ginv(t(x) %*%x))*varhat[j])
    ttest[,j] <- BBhat[,j]/covbb[,j]
    pval[,j] <- 2*pt(abs(ttest[,j]),df=(nrowy-m*ncolx), lower = F)
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

rownames(BBhat) <- c(colnames(x), colnames(x))
colnames(BBhat) <- colnames(y)
rownames(bhat1) <- c(colnames(x), colnames(x))
colnames(bhat1) <- colnames(y)
results <- list(BBhat, covbb, ttest, pval, cgam1, omega[[iter]], residuals1, bhat1)
names(results) <- c('Bhat','St.Dev.', 't-test', 'pval', 'C-gamma', 'Omega', 'Residuals', 'Output')
  return(results)
}
