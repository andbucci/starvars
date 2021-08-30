VLSTAR <- function(y, exo = NULL, p = 1,
                   m = 2, st = NULL, constant = TRUE, starting = NULL,
                   method = c('ML', 'NLS'), n.iter = 500,
                   singlecgamma = FALSE,
                   epsilon = 10^(-3)){
  y <- as.matrix(y)
  x <- exo
  ncores <- detectCores()
  method <- match.arg(method)
  ##Checks and warnings
  if (anyNA(y))
    stop("\nNAs in y.\n")
  if(m < 2)
    stop('The number of regimes should be greater than one.')
  if(is.null(st))
    stop('The transition variable must be supplied.')
  if(is.null(x)){
    if(length(y[,1]) != length(st))
      stop('The length of the variables does not match!')
  }else{
    if(length(y[,1]) != length(as.matrix(x[,1])) | length(st) != length(as.matrix(x[,1])) | length(y[,1]) != length(st))
      stop('The length of the variables does not match!')
  }
  if(is.null(starting)){
    stop('Starting values should be provided.')
  }
  if(!is.list(starting)){
    stop('Starting c and gamma should be put in a list.')
  }
  if(is.null(p) || p < 1){
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
  ylag <- stats::lag(yt, -(1:p))
  ylag <- as.matrix(ylag)
  if(p>1){
    lagg <- p-1
    ylag <- ylag[-(1:lagg),]
  }
  y <- y[-c(1:p), ]
  ncoly <- ncol(y)
  if(!is.null(starting)){
    if(length(starting)!= (m-1)){
      stop('The length of the list of initial values should be equal to m-1.')
    }else{
      if(any(unlist(lapply(starting, ncol))!=2) | any(unlist(lapply(starting, nrow))!=ncoly)){
        stop('Each element of the starting argument should have two columns and n rows.')
      }

    }
  }
  ncolylag <- ncoly*p
  nrowy <- nrow(y)
  ncolx1 <- ncol(x)
  const <- rep(1, nrowy)
  if (constant == TRUE){
    if(!is.null(exo)){
      x1a <- as.matrix(x[-c(1:p),])
      x <- as.matrix(cbind(const,ylag,x1a))
    }else{
      x <- as.matrix(cbind(const,ylag))
    }
    ncolx <- ncol(x)
  }  else{
    x <- as.matrix(x)
  }
  nrowx <- nrow(x)
  st <- st[(1+p):length(st)]
  ncolx <- ncol(x)
  param.init <- list()
  if(singlecgamma == TRUE){
    param.init$gamma <- 1
    param.init$c <- mean(y)
  } else{
    param.init$gamma <- rep(1L, ncoly)
    param.init$c <- colMeans(y)
  }
  ny <- ifelse(singlecgamma == TRUE, 1, ncoly)
  q <- ncol(x)-ncolylag
  PARAM <- starting

####Estimating VLSTAR model####

  ##Estimating initial values to be used in the iterative algorithm
  In <- diag(ncoly)
  glog <- matrix(ncol=ny, nrow = nrowy)
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
    for (t in 1:(m-1)){
      for (j in 1 : ny){
        glog[i,j] <- (1L+exp(-PARAM[[t]][j,1]*(st[i]-PARAM[[t]][j,2])))^(-1)
      }
      if(singlecgamma == TRUE){
        Gt <- diag(rep(glog[i,1], ncoly))
      }else{
        Gt <- diag(glog[i,])
      }
      GT[[t]] <- Gt
    }
    Gtilde[[i]] <- t(cbind(In, do.call(cbind,GT)))
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
  BB <- invvec(Bhat, ncol = (m*ncoly), nrow = (ncolylag + q)) ##Estimated coefficients
  ##Calculating the estimated covariance matrix
  resi <- list()
  Ehat1 <- matrix(NA, ncol = ncoly, nrow = nrowy)
  for (i in 1:nrowx){
    resi[[i]] <- y[i, ] - t(Gtilde[[i]])%*%t(BB)%*%x[i,]
  }
  Ehat1 <- t(do.call("cbind", resi))
  #Estimated covariance matrix
  Omegahat <- (t(Ehat1)%*%Ehat1)/(nrowy)

  data = list(y = y, x = x, m = m, BB = BB, Omegahat = Omegahat, st = st, singlecgamma = singlecgamma)

  #Inizialization of iter
  iter <- 0
  ll0 <- 10^(4)
  epsi <- epsilon ##Value used as convergence check
  loglik1 <- NULL
  bbhat <- list()
  omega <- list()
  PARAM1 <- list()
  for(t in 1:(m-1)){
    PARAM1[[t]] <- as.data.frame(PARAM[[t]])
  }
  param <- rbindlist(PARAM1)
  param <- vec(as.matrix(param)) ##the parameters c and gamma are vectorized in order to be passed in the optimParallel function
  err <- 10^5

  #Log-likelihood to be optimized in the ML method and used to check convergence in both methods
  loglike <- function(param, data){
    m = data$m
    y = data$y
    x = data$x
    st = data$st
    BB = data$BB
    Omegahat = data$Omegahat
    ncoly = ncol(y)
    nrowy = nrow(y)
    singlecgamma = data$singlecgamma
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
        if(singlecgamma == TRUE){
          Gt <- diag(rep(glog[z,1], ncoly))
        }else{
          Gt <- diag(glog[z,])
        }
        GT[[t]] <- Gt
      }
      Gtilde[[z]] <- t(cbind(In, do.call(cbind,GT)))
      dify[z] <-  t(y[z, ] - t(Gtilde[[z]])%*%t(BB)%*%x[z,])%*%MASS::ginv(Omegahat)%*%(y[z, ] - t(Gtilde[[z]])%*%t(BB)%*%x[z,])
    }
    sumdif <- sum(dify)
    logll <- -(nrowy*log(det(Omegahat))/2L) - sumdif/2L  - (nrowy*ncoly/2L)*log(2L*pi)##Normal distribution assumed
    return(-logll)
  }

loglike2 <- function(y, resid1, omega){
    nrowy <- nrow(as.matrix(y))
    logll <- -(nrowy/2)*log(2*pi) -(nrowy/2)*log(omega) - (t(resid1)%*%resid1)/(2*omega)
    return(logll)
}
##Actual iteration to estimate the coefficients
if(method == 'ML'){
    #NLS Estimation of Bhat and Omegahat to be used in the first iteration of maximum likelihood
  message('Maximum likelihood estimation\n')
  errdif <- 10^5
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
      cl <- makeCluster(ncores)     # set the number of processor cores
      setDefaultCluster(cl=cl)
      param1 <- optimParallel(par = as.vector(param), fn = loglike, lower = c(low1, apply(y, 2, min)),
                      method="L-BFGS-B", data = data)

      cgam1 <- matrix(param1$par, ncol = 2, nrow = (ny*(m-1)))

      #2.Maximum likelihood estimation of Bhat with new values of gamma and c
      glog <- matrix(ncol=ny, nrow = nrowy)
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
          for (j in 1 : ny){
            glog[i,j] <- (1+exp(-cgam1[j,1]*(st[i]-cgam1[j,2])))^(-1)
          }
          if(singlecgamma == TRUE){
            Gt <- diag(rep(glog[i,1], ncoly))
          }else{
            Gt <- diag(glog[i,])
          }
          GT[[t]] <- Gt
        }
        Gtilde[[i]] <- t(cbind(In, do.call(cbind,GT)))
        GG[[i]] <- Gtilde[[i]]%*%t(Gtilde[[i]])
        XX[[i]] <- x[i,] %*%t(x[i,])
        XY[[i]] <- x[i, ]%*%t(y[i,])
        XYOG[[i]] <- vec(XY[[i]]%*%ginv(Omegahat)%*%t(Gtilde[[i]]))
        PsiOmegaPsi[[i]] <- Gtilde[[i]]%*%ginv(Omegahat)%*%t(Gtilde[[i]])
        kro[[i]] <- kronecker(PsiOmegaPsi[[i]], XX[[i]])
      }
      xyog <- Reduce(`+`, XYOG)/nrowy
      kroxx <- Reduce(`+`, kro)/nrowy
      Bhat <- t(t(xyog)%*%kroxx)
      BB <- invvec(Bhat, ncol = (m*ncoly), nrow = (ncolylag + q))
      resi <- list()
      fitte <- matrix(nrow = nrowy, ncol = ncoly)
      Ehat1 <- matrix(NA, ncol = ncoly, nrow = nrowy)
      for (i in 1:nrowx){
        resi[[i]] <- y[i, ] - t(Gtilde[[i]])%*%t(BB)%*%x[i,]
        fitte[i,] <- t(t(Gtilde[[i]])%*%t(BB)%*%x[i,])
      }
      Ehat1 <- t(do.call("cbind", resi))
      Omegahat <- (t(Ehat1)%*%Ehat1)/(nrowy)
      data$Omegahat = Omegahat
      data$BB = BB
      #3. Convergence check
      ll1 <- loglike(param = param1$par, data)
      if((ll1-ll0)>0){
        param <- as.matrix(vec(cgam1))
      }
      err <- abs(ll1 - ll0)
      ll0 <- ll1
      loglik1[iter] <- ll1
      bbhat[[iter]] <- BB
      omega[[iter]] <- Omegahat

      message(paste("iteration", iter, "complete\n"))

      cat(paste('Log-likelihood:', round(ll0,3), '\n'))

      if (err<epsi | iter == n.iter) message('Converged\n')}
  } else{
    #NLS Estimation of Bhat and Omegahat to be used in the first iteration of minimizing Qt

    message('NLS estimation\n')

    #Inizialization of iter
    #Sum of squared error to be optimized
    ssq1 <- function(param, data){
      y <- data$y
      ncoly = ncol(y)
      nrowy = nrow(y)
      m = data$m
      BB = data$BB
      x = data$x
      st = data$st
      singlecgamma = data$singlecgamma
      gamma <- param[1:((m-1)*ncoly)]
      c <- param[(ncoly*(m-1)+1):length(param)]
      gamma1 <- matrix(gamma, ncol = (m-1))
      c1 <- matrix(c, ncol = (m-1))
      glog <- matrix(ncol=ncoly, nrow = nrowy)
      GT <- list()
      Gtilde <- list()
      dify <- matrix(ncol = 1, nrow = nrowy)
      ncoly <- dim(y)[2]
      In <- diag(ncoly)
      for (z in 1:nrow(y)){
        for(t in 1:(m-1)){
          for (o in 1:ncoly){
            gammao <- gamma1[o,t]
            co <- c1[o,t]
            glog[z,o] <- (1L+exp(-gammao*(st[z]-co)))^(-1)}
          if(singlecgamma == TRUE){
            Gt <- diag(rep(glog[z,1], ncoly))
          }else{
            Gt <- diag(glog[z,])
          }
          GT[[t]] <- Gt
        }
        Gtilde[[z]] <- t(cbind(In, do.call(cbind,GT)))
        dify[z] <-  t(y[z, ] - t(Gtilde[[z]])%*%t(BB)%*%x[z,])%*%(y[z, ] - t(Gtilde[[z]])%*%t(BB)%*%x[z,])
      }
      sumdif <- sum(dify)
      return(sumdif)
    }
    #Convergence algorithm
    #1. NLS estimation of gamma and c with NLS estimates of Bhat and Omegahat
    #2. NLS of Bhat with new values of gamma and c
    #3. Convergence check
    #4. 1-2-3 until convergence
    while (iter < n.iter & err > epsi){
      Sys.sleep(0)
      iter <- iter+1L
      #Parameters
      low1 <- replicate(ny, 0)
      #1.Maximum likelihood estimation of gamma and c
      cl <- makeCluster(ncores)     # set the number of processor cores
      setDefaultCluster(cl=cl)
      param1 <- optimParallel(par = as.vector(param), fn = ssq1, lower = c(low1, apply(y, 2, min)),
                      method="L-BFGS-B", data = data)

      cgam1 <- matrix(param1$par, ncol = 2L, nrow = (ny*(m-1)))

      #2.NLS estimation of Bhat with new values of gamma and c
      glog <- matrix(ncol=ny, nrow = nrowy)
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
          for (j in 1 : ny){
            glog[i,j] <- (1+exp(-cgam1[((t-1)*ny + j),1]*(st[i]-cgam1[((t-1)*ny + j),2])))^(-1)
          }
          if(singlecgamma == TRUE){
            Gt <- diag(rep(glog[i,1], ncoly))
          }else{
            Gt <- diag(glog[i,])
          }
          GT[[t]] <- Gt
        }
        Gtilde[[i]] <- t(cbind(In, do.call(cbind,GT)))
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
      BB <- invvec(Bhat, ncol = (m*ncoly), nrow = (ncolylag + q))
      resi <- list()
      resiresi <- list()
      fitte <- matrix(nrow = nrowy, ncol = ncoly)
      for (o in 1:nrowx){
        resi[[o]] <- y[o, ] - t(Gtilde[[o]])%*%t(BB)%*%x[o,]
        resiresi[[o]] <- resi[[o]]%*%t(resi[[o]])
        fitte[o,] <- t(t(Gtilde[[o]])%*%t(BB)%*%x[o,])
      }
      Ehat1 <- Reduce("+", resiresi)
      Omegahat <- Ehat1/(nrowy-1L)

      data$Omegahat = Omegahat
      data$BB = BB

      #3. Convergence check
      ll1 <- loglike(param = param1$par, data)
      if((ll1-ll0)>0){
        param <- as.matrix(vec(cgam1))
      }
      err <- abs(ll1 - ll0)
      ll0 <- ll1
      loglik1[iter] <- ll1
      bbhat[[iter]] <- BB
      omega[[iter]] <- Omegahat

      message(paste("iteration", iter, "complete\n"))

      cat(paste('Log-likelihood:',round(ll0, 3), '\n'))

      if(iter > 50){
        if(loglik1[[iter]] == loglik1[[(iter-2)]]){
          iter <- n.iter
        }}

      if (err<epsi | iter == n.iter) message('Converged\n')}
  }

    #Calculating residuals and estimating standard errors
    residuals1 <- t(do.call("cbind", resi))
    varhat <- diag(omega[[iter]])
    bb1 <- bbhat[[iter]][,1:ncoly]
    bb2 <- list()
    for(t in 1:(m-1)){
      bb2[[t]] <- as.data.frame(bbhat[[iter]][,(ncoly*(t-1)+1):(ncoly*t)] +
                                  bbhat[[iter]][,(ncoly*(t)+1):(ncoly*(t+1))])
    }
    bb4 <- as.matrix(rbindlist(bb2))
    BBhat <- rbind(bb1, bb4) ##Estimates of the coefficients

    ##Calculating standard errors, t-test and p-values
    colnames(cgam1) <- c('gamma', 'c')
    covbb <- matrix(nrow = m*ncolx, ncol = ncoly)
    ttest <- matrix(nrow = m*ncolx, ncol = ncoly)
    pval <- matrix(nrow = m*ncolx, ncol = ncoly)
    ee <- matrix(nrow = m*ncolx, ncol = ncoly)
    for (j in 1 : ncoly){
      covbb[,j] <- sqrt(diag(ginv(t(x) %*%x))*varhat[j])
      ttest[,j] <- BBhat[,j]/covbb[,j]
      pval[,j] <- 2*pt(abs(ttest[,j]),df=(nrowy*m-m*(ncolx)), lower.tail = FALSE)
    }
    ##Significances used for the summary function
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

    ##Calculating univariate log-likelihoods, AIC and BIC criteria
    residui <- residuals1
    omega1 <- diag(omega[[iter]])
    k <- nrow(BBhat)
    ll2 <- NULL
    AIC1 <- NULL
    BIC1 <- NULL
    for (l in 1:ncoly){
      ll2[l] <- loglike2(y[,l], residui[,l], omega1[l])
      AIC1[l] <- 2*k - 2*ll2[l]
      BIC1[l] <- -2*ll2[l] + k*log(nrowy)
    }
    names1 <- list()
    for(j in 1:m){
      names1[[j]] <- as.data.frame(paste(colnames(x), 'm_', j))
    }
    names1 <- as.matrix(rbindlist(names1))
    rownames(BBhat) <- names1
    colnames(BBhat) <- colnames(y)
    modeldata <- list(y, x)
    fitte <- fitte[!is.na(fitte[,1]),]
    results <- list(BBhat, covbb, ttest, pval, cgam1, omega[[iter]], fitte, residuals1, ll1, ll2, AIC1, BIC1, Gtilde, modeldata, BB, m, p,
                    st, y, exo, constant, method, singlecgamma)
    names(results) <- c('Bhat','StDev', 'ttest', 'pval', 'Gammac', 'Omega', 'fitted', 'residuals', 'MultiLL', 'LL', 'AIC',
                        'BIC', 'Gtilde', 'Data', 'B', 'm', 'p', 'st', 'yoriginal', 'exo', 'constant', 'method', 'singlecgamma')
  class(results) = 'VLSTAR'
  return(results)
}


