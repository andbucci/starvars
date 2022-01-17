#' VLSTAR Prediction
#'
#' One-step or multi-step ahead forecasts, with interval forecast, of a VLSTAR object.
#'
#' @param object An object of class \sQuote{\code{VLSTAR}} obtained through \command{VLSTAR()}
#' @param n.ahead An integer specifying the number of ahead predictions
#' @param conf.lev Confidence level of the interval forecast
#' @param st.new Vector of new data for the transition variable
#' @param M An integer with the number of errors sampled for the \code{Monte Carlo} method
#' @param B An integer with the number of errors sampled for the \code{bootstrap} method
#' @param st.num An integer with the index of dependent variable if \code{st.new} is \code{NULL}
#' and the transition variable is a lag of one of the dependent variables
#' @param method A character identifying which multi-step ahead method should be used among \code{naive}, \code{Monte Carlo} and \code{bootstrap}
#' @param newdata \code{data.frame} or \code{matrix} of new data for the exogenous variables
#' @param \dots further arguments to be passed to and from other methods
#' @return A \code{list} containing:
#'\item{forecasts}{\code{data.frame} of predictions for each dependent variable and the (1-\eqn{\alpha}) prediction intervals}
#'\item{y}{a matrix of values for y}
#' @aliases predict print.vlstarpred
#' @references Granger C.W.J. and Terasvirta T. (1993), Modelling Non-Linear Economic Relationships. \emph{Oxford University Press};
#'
#' Lundbergh S. and Terasvirta T. (2007), Forecasting with Smooth Transition Autoregressive Models. \emph{John Wiley and Sons};
#'
#' Terasvirta T. and Yang Y. (2014), Specification, Estimation and Evaluation of Vector Smooth Transition Autoregressive Models with Applications. \emph{CREATES Research Paper 2014-8}
#' @author Andrea Bucci and Eduardo Rossi
#' @seealso \code{\link{VLSTAR}} for log-likehood and nonlinear least squares estimation of the VLSTAR model.
#' @keywords VLSTAR



predict.VLSTAR <- function(object, ..., n.ahead = 1, conf.lev = 0.95, st.new = NULL, M = 5000, B = 1000,
                           st.num = NULL, newdata = NULL, method = c('naive', 'Monte Carlo', 'bootstrap')){
  st_new <- st.new
  method <- match.arg(method)
  resid <- object$residuals
  ncoly <- ncol(object$Data[[1]])
  nrowy <- nrow(object$Data[[1]])
  alpha = 1-conf.lev
  k <- ncol(object$exo)
  In <- diag(ncol(object$Data[[1]]))
  datamat <- object$Data[[2]]
if(is.null(object$exo)){
  if(object$constant == T){
    constant <- matrix(rep(1, n.ahead), nrow = n.ahead, ncol = 1)
    yy <- datamat[nrow(datamat),c(2:(ncoly*object$p+1))]
  }else{
    yy <- datamat[nrow(datamat),c(1:(ncoly*object$p))]
  }
  }else{
    if(is.null(newdata) | length(newdata) != k){
      stop('Please, provide valid new data for the exogenous variables')
    }
  if(object$constant == TRUE){
   constant <- matrix(rep(1, n.ahead), nrow = n.ahead, ncol = 1)
   yy <- datamat[nrow(datamat),c(2:(ncoly*object$p+1))]
    }else{
      yy <- datamat[nrow(datamat),c(1:(ncoly*object$p))]
  }

newdata <- as.matrix(newdata)
}

  BB <- object$B
  BB1 <- object$B[,seq(1, ncol(BB))[c(rep(TRUE, ncoly), rep(FALSE, ncoly))]]
  BB2 <- object$B[,seq(1, ncol(BB))[c(rep(FALSE, ncoly), rep(TRUE, ncoly))]]
  Gtilde <- t(object$Gtilde[[length(object$Gtilde)]])
  pred <- matrix(NA, ncol = ncoly, nrow = n.ahead)
  Z <- c(constant[1, ], yy, newdata)
  pred[1,] <- as.matrix((Gtilde%*%t(BB)%*%Z))
  tmp <- pred[1, ]
  yy <- c(tmp)
  if(is.null(st.new)){
    if(is.null(st.num)){
      stop('Please, provide valid index for the transition variable')
    }

    st_new <-c(as.matrix(object$st[length(object$st)]), pred[1,st.num])
  }
if(!is.null(st.new)){
  if(length(st.new) != n.ahead){
    stop('The length of the new data for the transition variable should be equal to the number of steps ahead')
}
}

  lower <- matrix(NA, ncol = ncoly, nrow = n.ahead)
  upper <- matrix(NA, ncol = ncoly, nrow = n.ahead)
  for(j in 1:ncoly){
    std_err <- sqrt((ncoly-object$p-1)/(nrowy-object$p-ncoly)*qf(1-alpha/2, df1 = ncoly, df2 = (nrowy-object$p-ncoly)))*sqrt((1+t(Z)%*%MASS::ginv(t(object$Data[[2]])%*%object$Data[[2]])%*%Z)*object$Omega[j,j])
    lower[1,j] <- pred[1,j]-std_err
    upper[1,j] <- pred[1,j]+std_err
  }

  if(n.ahead >1){
  if(method == 'naive'){
    for(i in 2 : n.ahead){
      Z <- c(constant[i, ], yy, newdata)
      GT <- matrix(NA, ncol = 1, nrow = ncoly)
      for(j in 1:ncoly){
        GT[j,] <- (1+exp(-object$Gammac[j,1]*(st_new[i]-object$Gammac[j,2])))^(-1)
      }
      Gt <- diag(GT[,1])
      Gtilde <- t(cbind(In, Gt))
      pred[i,] <- as.matrix((t(Gtilde)%*%t(BB)%*%Z))
      tmp <- pred[i, ]
      yy <- c(tmp)
      if(is.null(st.new)){
        st_new <- c(as.vector(st_new),tmp[st.num])
      }
      for(j in 1:ncoly){
        std_err <- sqrt((ncoly-object$p-1)/(nrowy-object$p-ncoly)*qf(1-alpha/2, df1 = ncoly, df2 = (nrowy-object$p-ncoly)))*sqrt((1+t(Z)%*%MASS::ginv(t(object$Data[[2]])%*%object$Data[[2]])%*%Z)*object$Omega[j,j])
        lower[i,j] <- pred[1,j]-std_err
        upper[i,j] <- pred[1,j]+std_err
      }
    }

  }else if(method == 'Monte Carlo'){
    for(i in 2 : n.ahead){
    epsilon_mc <- matrix(NA, nrow = M, ncol = ncoly)
    for(k in 1:ncoly){
      epsilon_mc[,k] <- rnorm(M,mean=0,sd=sd(resid[,k])*3)
    }
    nonlinear_MC <- matrix(NA, ncol = M, nrow = ncoly)
      for(m in 1 : M){
        ymc <- yy + epsilon_mc[m, ]
        Zmc <- c(constant[i, ], ymc, newdata)
        GT <- matrix(NA, ncol = 1, nrow = ncoly)
        for(j in 1:ncoly){
          GT[j,] <- (1+exp(-object$Gammac[j,1]*(st_new[i]-object$Gammac[j,2])))^(-1)
        }
        Gt <- diag(GT[,1])
        Gtilde <- t(cbind(In, Gt))
        nonlinear_MC[,m] <- t(Gtilde)%*%t(BB)%*%as.matrix(Zmc)
      }
    nonlinearMC <- t(nonlinear_MC)
      Z <- c(constant[i, ], yy, newdata)
      pred[i,] <- t(t(BB1)%*%Z + rowMeans(nonlinear_MC))
      tmp <- pred[i, ]
      yy <- c(tmp)
      if(is.null(st.new)){
        st_new <- c(as.vector(st_new),tmp[st.num])
      }
      for(j in 1:ncoly){
        lower[i,j] <- t(t(BB1)%*%Z)[j] + quantile(nonlinearMC[,j], (alpha/2))
        upper[i,j] <-t(t(BB1)%*%Z)[j] + quantile(nonlinearMC[,j], (1-alpha/2))
      }
    }
  }else if(method == 'bootstrap'){
    for(i in 2 : n.ahead){
      epsilon_bo <- matrix(NA, nrow = B, ncol = ncoly)
      nonlinear_BO <- matrix(NA, ncol = B, nrow = ncoly)
      for(k in 1:ncoly){
        for(b in 1 : B){
        epsilon_bo[b,k] <- sample(resid[,k],1)
        ybo <- yy + epsilon_bo[b, ]
        Zbo <- c(constant[i, ], ybo, newdata)
        GT <- matrix(NA, ncol = 1, nrow = ncoly)
        for(j in 1:ncoly){
          GT[j,] <- (1+exp(-object$Gammac[j,1]*(st_new[i]-object$Gammac[j,2])))^(-1)
        }
        Gt <- diag(GT[,1])
        Gtilde <- t(cbind(In, Gt))
        nonlinear_BO[,b] <- t(Gtilde)%*%t(BB)%*%Zbo
      }}
      nonlinearBO <- t(nonlinear_BO)
      Z <- c(constant[i, ], yy, newdata)
      pred[i,] <- t(t(BB1)%*%Z + rowMeans(nonlinear_BO))
      tmp <- pred[i, ]
      yy <- c(tmp)
      if(is.null(st.new)){
        st_new <- c(as.vector(st_new),tmp[st.num])
      }
      for(j in 1:ncoly){
        lower[i,j] <- t(t(BB1)%*%Z)[j] + quantile(nonlinearBO[,j], (alpha/2))
        upper[i,j] <-t(t(BB1)%*%Z)[j] + quantile(nonlinearBO[,j], (1-alpha/2))
      }
    }
    }
}

colnames(pred) <- colnames(object$Bhat)
rownames(pred) <- paste("Step", 1:n.ahead)


forecasts <- list()
for(i in 1 : ncoly){
  forecasts[[i]] <- cbind(pred[, i], lower[, i], upper[, i])
  colnames(forecasts[[i]]) <- c("fcst", paste("lower ", (conf.lev*100), "%", sep = ''), paste("upper ", (conf.lev*100), "%", sep = ''))
}
names(forecasts) <- colnames(object$Bhat)
results <- list(forecasts = forecasts, y = object$Data[[1]])
class(results) <- "vlstarpred"
  return(results)
}
