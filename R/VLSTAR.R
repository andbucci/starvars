#' VLSTAR- Estimation
#'
#' This function allows the user to estimate the coefficients of a VLSTAR model with \emph{m} regimes through maximum likelihood or nonlinear least squares.
#' The set of starting values of Gamma and C for the convergence algorithm can be either passed or obtained via searching grid.
#'
#' The multivariate smooth transition model is an extension of the smooth transition regression model introduced by Bacon and Watts (1971)  (see also Anderson and Vahid, 1998). The general model is
#' \deqn{y_{t} = \mu_0+\sum_{j=1}^{p}\Phi_{0,j}\,y_{t-j}+A_0 x_t \cdot G_t(s_t;\gamma,c)[\mu_{1}+\sum_{j=1}^{p}\Phi_{1,j}\,y_{t-j}+A_1x_t]+\varepsilon_t}
#' where \eqn{\mu_{0}} and \eqn{\mu_{1}} are the \eqn{\tilde{n} \times 1} vectors of intercepts, \eqn{\Phi_{0,j}} and \eqn{\Phi_{1,j}} are square
#' \eqn{\tilde{n}\times\tilde{n}} matrices of parameters for lags \eqn{j=1,2,\dots,p}, A_0 and A_1 are \eqn{\tilde{n}\times k} matrices of parameters,
#' x_t is the \eqn{k \times 1} vector of exogenous variables and \eqn{\varepsilon_t} is the innovation. Finally, \eqn{G_t(s_t;\gamma,c)} is a \eqn{\tilde{n}\times \tilde{n}} diagonal matrix of transition function at time \emph{t}, such that
#' \deqn{G_t(s_t;\gamma,c)=\{G_{1,t}(s_{1,t};\gamma_{1},c_{1}),G_{2,t}(s_{2,t};\gamma_{2},c_{2}),
#' \dots,G_{\tilde{n},t}(s_{\tilde{n},t};\gamma_{\tilde{n}},c_{\tilde{n}})\}.}
#'
#' Each diagonal element \eqn{G_{i,t}^r} is specified as a logistic cumulative density functions, i.e.
#' \deqn{G_{i,t}^r(s_{i,t}^r; \gamma_i^r, c_i^r) = \left[1 + \exp\big\{-\gamma_i^r(s_{i,t}^r-c_i^r)\big\}\right]^{-1}}
#' for \eqn{latex}{i = 1,2, \dots, \tilde{n}} and \eqn{r=0,1,\dots,m-1}, so that the first model is a Vector Logistic Smooth Transition AutoRegressive
#' (VLSTAR) model.
#' The ML estimator of \eqn{\theta} is obtained by solving the optimization problem
#' \deqn{\hat{\theta}_{ML} = arg \max_{\theta}log L(\theta)}
#' where \eqn{log L(\theta)} is the log-likelihood function of VLSTAR model, given by
#' \deqn{  ll(y_t|I_t;\theta)=-\frac{T\tilde{n}}{2}\ln(2\pi)-\frac{T}{2}\ln|\Omega|-\frac{1}{2}\sum_{t=1}^{T}(y_t-\tilde{G}_tB\,z_t)'\Omega^{-1}(y_t-\tilde{G}_tB\,z_t)}
#' The NLS estimators of the VLSTAR model are obtained by solving the optimization problem
#' \deqn{\hat{\theta}_{NLS} = arg \min_{\theta}\sum_{t=1}^{T}(y_t - \Psi_t'B'x_t)'(y_t - \Psi_t'B'x_t).}
#' Generally, the optimization algorithm may converge to some local minimum. For this reason, providing valid starting values of \eqn{\theta} is crucial. If there is no clear indication on the initial set of parameters, \eqn{\theta}, this can be done by implementing a grid search. Thus, a discrete grid in the parameter space of \eqn{\Gamma} and C is create to obtain the estimates of B conditionally on each point in the grid. The initial pair of \eqn{\Gamma} and C producing the smallest sum of squared residuals is chosen as initial values, then the model is linear in parameters.
#' The algorithm is the following:
#' \enumerate{
#' \item Construction of the grid for \eqn{\Gamma} and C, computing \eqn{\Psi} for each poin in the grid
#' \item Estimation of \eqn{\hat{B}} in each equation, calculating the residual sum of squares, \eqn{Q_t}
#' \item Finding the pair of \eqn{\Gamma} and C providing the smallest \eqn{Q_t}
#' \item Once obtained the starting-values, estimation of parameters, \emph{B}, via nonlinear least squares (NLS)
#' \item Estimation of \eqn{\Gamma} and C given the parameters found in step 4
#' \item Repeat step 4 and 5 until convergence.
#' }
#'
#' @param y \code{data.frame} or \code{matrix} of dependent variables of dimension \code{(Txn)}
#' @param exo (optional) \code{data.frame} or \code{matrix} of exogenous variables of dimension \code{(Txk)}
#' @param p lag order
#' @param m number of regimes
#' @param st single transition variable for all the equation of dimension \code{(Tx1)}
#' @param constant \code{TRUE} or \code{FALSE} to include or not the constant
#' @param starting set of intial values for Gamma and C, inserted as a list of length \code{m-1}.
#' Each element of the list should contain a \code{data.frame} with 2 columns (one for Gamma and one for c), and \code{n} rows.
#' @param method Fitting method: maximum likelihood or nonlinear least squares.
#' @param singlecgamma \code{TRUE} or \code{FALSE} to use single gamma and c
#' @param n.iter number of iteration of the algorithm until forced convergence
#' @param epsilon convergence check measure
#' @param ncores Number of cores used for parallel computation. Set to \code{NULL} by default and automatically calculated.
#' @return An object of class \code{VLSTAR}, with standard methods.
#' @references Anderson H.M. and Vahid F. (1998), Testing multiple equation systems for common nonlinear components. \emph{Journal of Econometrics}. 84: 1-36
#'
#' Bacon D.W. and Watts D.G. (1971), Estimating the transition between two intersecting straight lines. \emph{Biometrika}. 58: 525-534
#'
#' Terasvirta T. and Yang Y. (2014), Specification, Estimation and Evaluation of Vector Smooth Transition Autoregressive Models with Applications. \emph{CREATES Research Paper 2014-8}
#' @author Andrea Bucci
#' @export
#' @exportClass VLSTAR
#' @importFrom MASS ginv
#' @importFrom dplyr if_else
#' @importFrom ks invvec
#' @importFrom matrixcalc vec vech
#' @importFrom optimParallel optimParallel
#' @importFrom zoo zoo
#' @importFrom lessR to
#' @importFrom data.table rbindlist
#' @importFrom parallel detectCores makeCluster setDefaultCluster stopCluster
#' @importFrom stats lag residuals
#' @keywords VLSTAR
#' @examples
#' \donttest{
#' data(Realized)
#' y <- Realized[-1,1:10]
#' y <- y[-nrow(y),]
#' st <- Realized[-nrow(Realized),1]
#' st <- st[-length(st)]
#' stvalues <- startingVLSTAR(y, p = 1, n.combi = 3,
#'  singlecgamma = FALSE, st = st)
#'fit.VLSTAR <- VLSTAR(y, p = 1, singlecgamma = FALSE, starting = stvalues,
#'  n.iter = 3, st = st, method ='NLS')
#'# a few methods for VLSTAR
#'print(fit.VLSTAR)
#'summary(fit.VLSTAR)
#'plot(fit.VLSTAR)
#'predict(fit.VLSTAR, st.num = 1, n.ahead = 1)
#'logLik(fit.VLSTAR, type = 'Univariate')
#'coef(fit.VLSTAR)}
#'
VLSTAR <- function(y, exo = NULL, p = 1,
                   m = 2, st = NULL, constant = TRUE, starting = NULL,
                   method = c('ML', 'NLS'), n.iter = 500,
                   singlecgamma = FALSE,
                   epsilon = 10^(-3), ncores = NULL){
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if(is.null(ncores)){
    if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    ncores <- 2L
  } else {
    ncores <- detectCores()
  }}
  y <- as.matrix(y)
  x <- exo
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
  kro <- list()
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
    kro[[i]] <- kronecker(Gtilde[[i]], x[i,])
  }
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
  omega <- list()
  PARAM1 <- list()
  for(t in 1:(m-1)){
    PARAM1[[t]] <- as.data.frame(PARAM[[t]])
  }
  param <- rbindlist(PARAM1)
  param <- vec(as.matrix(param)) ##the parameters c and gamma are vectorized in order to be passed in the optimParallel function
  err <- 10^5


##Actual iteration to estimate the coefficients
cl <- makeCluster(ncores)     # set the number of processor cores
setDefaultCluster(cl=cl)
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
      param1 <- optimParallel(par = as.vector(param), fn = loglike, lower = c(low1, apply(y, 2, min)),
                      data = data, parallel = list(cl = cl, forward = FALSE, loginfo = FALSE))
      cgam1 <- matrix(param1$par, ncol = 2, nrow = (ny*(m-1)))

      #2.Maximum likelihood estimation of Bhat with new values of gamma and c
      glog <- rep(0, ncol(y))
      GT <- list()
      Gtilde <- list()
      XX <- list()
      XYOG <- list()
      kro <- list()
      PsiOmegaPsi <- list()
      for (i in 1:nrowx){
        for(t in 1:(m-1)){
          for (j in 1 : ny){
            glog[j] <- (1+exp(-cgam1[j,1]*(st[i]-cgam1[j,2])))^(-1)
          }
          if(singlecgamma == TRUE){
            GT[[t]] <- diag(rep(glog[1], ncoly))
          }else{
            GT[[t]] <- diag(glog)
          }
        }
        Gtilde[[i]] <- t(cbind(In, do.call(cbind,GT)))
        XX[[i]] <- x[i,] %*%t(x[i,])
        XYOG[[i]] <- vec((x[i, ]%*%t(y[i,]))%*%ginv(Omegahat)%*%t(Gtilde[[i]]))
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

      message(paste("iteration", iter, "complete\n"))

      cat(paste('Log-likelihood:', round(ll0,3), '\n'))

      if (err<epsi | iter == n.iter) message('Converged\n')}
  } else{
    message('NLS estimation\n')

    #Inizialization of iter
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
      param1 <- optimParallel(par = as.vector(param), fn = SSQ, lower = c(low1, apply(y, 2, min)),
                      data = data, parallel = list(cl = cl, forward = FALSE, loginfo = FALSE))
      cgam1 <- matrix(param1$par, ncol = 2L, nrow = (ny*(m-1)))

      #2.NLS estimation of Bhat with new values of gamma and c
      glog <- rep(0, ncol(y))
      GT <- list()
      Gtilde <- list()
      kro <- list()
      for (i in 1:nrowx){
        for(t in 1:(m-1)){
          for (j in 1 : ny){
            glog[j] <- (1+exp(-cgam1[((t-1)*ny + j),1]*(st[i]-cgam1[((t-1)*ny + j),2])))^(-1)
          }
          if(singlecgamma == TRUE){
            GT[[t]] <- diag(rep(glog[1], ncoly))
          }else{
            GT[[t]] <- diag(glog)
          }
        }
        Gtilde[[i]] <- t(cbind(In, do.call(cbind,GT)))
        kro[[i]] <- kronecker(Gtilde[[i]], x[i,])
      }
      Bhat <- ginv(t(t(do.call("cbind", kro)))%*%(t(do.call("cbind", kro))))%*%t(t(do.call("cbind", kro)))%*%(vec(t(y)))
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

      message(paste("iteration", iter, "complete\n"))

      cat(paste('Log-likelihood:',round(ll0, 3), '\n'))

      if(iter > 50){
        if(loglik1[[iter]] == loglik1[[(iter-2)]]){
          iter <- n.iter
        }}

      if (err<epsi | iter == n.iter) message('Converged\n')}
  }
stopCluster(cl)

    #Calculating residuals and estimating standard errors
    residuals1 <- t(do.call("cbind", resi))
    varhat <- diag(Omegahat)
    bb1 <- BB[,1:ncoly]
    bb2 <- list()
    for(t in 1:(m-1)){
      bb2[[t]] <- as.data.frame(BB[,(ncoly*(t-1)+1):(ncoly*t)] +
                                  BB[,(ncoly*(t)+1):(ncoly*(t+1))])
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
    omega1 <- diag(Omegahat)
    k <- nrow(BBhat)
    ll2 <- NULL
    AIC1 <- NULL
    BIC1 <- NULL
    for (l in 1:ncoly){
      ll2[l] <- -(nrow(y)/2)*log(omega1[l]) - (t(residui[,l])%*%residui[,l])/(2*omega1[l])
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
    results <- list(BBhat, covbb, ttest, pval, cgam1, Omegahat, fitte, residuals1, ll1, ll2, AIC1, BIC1, Gtilde, modeldata, BB, m, p,
                    st, y, exo, constant, method, singlecgamma)
    names(results) <- c('Bhat','StDev', 'ttest', 'pval', 'Gammac', 'Omega', 'fitted', 'residuals', 'MultiLL', 'LL', 'AIC',
                        'BIC', 'Gtilde', 'Data', 'B', 'm', 'p', 'st', 'yoriginal', 'exo', 'constant', 'method', 'singlecgamma')
  class(results) = 'VLSTAR'
  return(results)
}


