#' Starting parameters for a VLSTAR model
#'
#' This function allows the user to obtain the set of starting values of Gamma and C for the convergence algorithm via searching grid.
#'
#' The searching grid algorithm allows for the optimal choice of the parameters \eqn{\gamma} and c by minimizing the sum of the Squared residuals for each possible combination.
#'
#' The parameter c is initialized by using the mean of the dependent(s) variable, while \eqn{\gamma} is sampled between 0 and 100.
#'
#' @param y \code{data.frame} or \code{matrix} of dependent variables of dimension \code{(Txn)}
#' @param exo (optional) \code{data.frame} or \code{matrix} of exogenous variables of dimension \code{(Txk)}
#' @param p lag order
#' @param m number of regimes
#' @param st single transition variable for all the equation of dimension \code{(Tx1)}
#' @param constant \code{TRUE} or \code{FALSE} to include or not the constant
#' @param n.combi Number of combination for the searching grid of Gamma and C
#' @param ncores  Number of cores used for parallel computation. Set to 2 by default
#' @param singlecgamma \code{TRUE} or \code{FALSE} to use single gamma and c
#' @return An object of class \code{startingVLSTAR}.
#' @importFrom foreach %dopar% foreach
#' @importFrom parallel setDefaultCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils setTxtProgressBar tail.matrix txtProgressBar
#' @seealso \code{\link{VLSTAR}}
#' @references Anderson H.M. and Vahid F. (1998), Testing multiple equation systems for common nonlinear components. \emph{Journal of Econometrics}. 84: 1-36
#'
#' Bacon D.W. and Watts D.G. (1971), Estimating the transition between two intersecting straight lines. \emph{Biometrika}. 58: 525-534
#'
#' Terasvirta T. and Yang Y. (2014), Specification, Estimation and Evaluation of Vector Smooth Transition Autoregressive Models with Applications. \emph{CREATES Research Paper 2014-8}
#' @author Andrea Bucci
#' @keywords VLSTAR
#' @export
#' @examples
#' \donttest{
#' data(Realized)
#' y <- Realized[-1,1:10]
#' y <- y[-nrow(y),]
#' st <- Realized[-nrow(Realized),1]
#' st <- st[-length(st)]
#' starting <- startingVLSTAR(y, p = 1, n.combi = 3,
#'                            singlecgamma = FALSE, st = st,
#'                            ncores = 1)}


startingVLSTAR <- function(y, exo = NULL, p = 1,
                   m = 2, st = NULL, constant = TRUE,
                   n.combi = NULL, ncores = 2,
                   singlecgamma = FALSE){
  y <- as.matrix(y)
  x <- exo
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
  if(is.null(n.combi)){
    stop('A number of combinations should be provided.')
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
  ylag <- as.matrix(stats::lag(yt, -(1:p)))
  if(p>1){
    lagg <- p-1
    ylag <- ylag[-(1:lagg),]
  }
  y <- y[-c(1:p), ]
  ncoly <- ncol(y)
  In <- diag(ncoly)
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
  COMBI <- list()
  GAMMA <- list()
  CJ <- list()
  for(t in 1:(m-1)){
      #Starting values for c and gamma
      gamma <- matrix(ncol = ny, nrow = n.combi)
      cj <- matrix(ncol = ny, nrow = n.combi)
      gamma[1,] <- param.init$gamma
      cj[1,] <- param.init$c
      GAMMA[[t]] <- gamma*(1+rnorm(1, mean = 0, sd = 5))
      CJ[[t]] <- cj*(1+rnorm(1, mean = 0, sd = 5))
      #Grid for c and gamma
      rangey <- matrix(nrow = ny, ncol = 2)
      combi <- list()
      for(j in 1:ny){
        rangey[j,] <- range(y[,j])
        CJ[[t]][2:n.combi,j] <- seq(from = rangey[j,1]*1.1, to = rangey[j,2], length.out = (n.combi-1))
        GAMMA[[t]][2:n.combi, j] <- seq(from = 0L, to = 100L, length.out = (n.combi-1))
        combi[[j]] <- expand.grid(CJ[[t]][,j], GAMMA[[t]][,j])
      }
      COMBI[[t]] <- combi
      }

#NLS for each combination of c and gamma

message(paste('Searching optimal c and gamma among', n.combi*n.combi, 'combinations\n'))
cl <- parallel::makeCluster(ncores)
doSNOW::registerDoSNOW(cl)
pb <- txtProgressBar(max = n.combi*n.combi, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
l = 1
ssq <- foreach(l = 1:(n.combi*n.combi), .errorhandling='pass', .combine = rbind,.options.snow = opts) %dopar% {
    Sys.sleep(0.1)
      glog <- matrix(ncol=ny, nrow = nrowy)
      GT <- list()
      Gtilde <- list()
      kro <- list()
      for (i in 1:nrowx){
        for(t in 1:(m-1)){
          for (j in 1 : ny){
            glog[i,j] <- (1+exp(-COMBI[[t]][[j]][l,2]*(st[i]-COMBI[[t]][[j]][l,1])))^(-1)
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
      Y <- matrixcalc::vec(t(y))
      Bhat <- MASS::ginv(t(M)%*%M)%*%t(M)%*%Y
      BB <- ks::invvec(Bhat, ncol = (m*ncoly), nrow = (ncolylag + q))
      resi <- list()
      Ehat <- matrix(NA, ncol = ncoly, nrow = nrowy)
      for (o in 1:nrowx){
        resi[[o]] <- y[o, ] - t(Gtilde[[o]])%*%t(BB)%*%x[o,]
      }
      Ehat <- t(do.call("cbind", resi))
      SSQ <- colSums(Ehat^2)
      return(SSQ)
}
close(pb)
stopCluster(cl)
#c and gamma minimazing the sum of squared residuals for each equation

cgamma <- matrix(nrow = ny, ncol = 2)
CGAMMA <- list()
for (t in 1:(m-1)){
  if(singlecgamma == T){
    cgamma[j,] <-  as.matrix(COMBI[[t]][[j]][which.min(rowSums(ssq)),])
    } else{
      for(j in 1:ny){
        cgamma[j,] <- as.matrix(COMBI[[t]][[j]][which.min(ssq[,j]),])
      }
      }
  CGAMMA[[t]] <- cgamma
  }


#Definition of c0 gamma0
PARAM <- list()
for (t in 1:(m-1)){
  cj <- as.matrix(CGAMMA[[t]][,1])
  gamma <- as.matrix(CGAMMA[[t]][,2])
  PARAM[[t]] <- cbind(gamma, cj)
  }
results <- PARAM
class(results) = 'startingVLSTAR'
return(results)
}


