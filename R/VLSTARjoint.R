#' Joint linearity test
#'
#' This function allows the user to test linearity against a Vector Smooth Transition Autoregressive Model with a single transition variable.
#'
#' Given a VLSTAR model with a unique transition variable, \eqn{s_{1t} = s_{2t} = \dots = s_{\widetilde{n}t} = s_t}, a generalization of the
#' linearity test presented in Luukkonen, Saikkonen and Terasvirta (1988) may be implemented.
#'
#' Assuming a 2-state VLSTAR model, such that
#' \deqn{y_t = B_{1}z_t + G_tB_{2}z_t + \varepsilon_t.}
#' Where the null \eqn{H_{0} : \gamma_{j} = 0}, \eqn{j = 1, \dots, \widetilde{n}}, is such that \eqn{G_t \equiv (1/2)/\widetilde{n}} and the
#' previous Equation is linear. When the null cannot be rejected, an identification problem of the parameter \eqn{c_{j}} in the transition function
#' emerges, that can be solved through a first-order Taylor expansion around \eqn{\gamma_{j} = 0}.
#'
#' The approximation of the logistic function with a first-order Taylor expansion is given by
#' \deqn{G(s_t; \gamma_{j},c_{j}) = (1/2) + (1/4)\gamma_{j}(s_t-c_{j}) + r_{jt}}
#' \deqn{= a_{j}s_t + b_{j} + r_{jt}}
#' where \eqn{a_{j} = \gamma_{j}/4}, \eqn{b_{j} = 1/2 - a_{j}c_{j}} and \eqn{r_{j}} is the error of the approximation. If \eqn{G_t} is specified as follows
#' \deqn{G_t = diag\big\{a_{1}s_t + b_{1} + r_{1t}, \dots, a_{\widetilde{n}}s_t+b_{\widetilde{n}} + r_{\widetilde{n}t}\big\}}
#' \deqn{= As_t + B + R_t}
#' where \eqn{A = diag(a_{1}, \dots, a_{\widetilde{n}})}, \eqn{B = diag(b_{1},\dots, b_{\widetilde{n}})} e \eqn{R_t = diag(r_{1t}, \dots,
#'                                                                                                                        r_{\widetilde{n}t})}, \eqn{y_t} can be written as
#'                                                                                                                        \deqn{y_t = B_{1}z_t + (As_t + B + R_t)B_{2}z_t+\varepsilon_t}
#'                                                                                                                        \deqn{= (B_{1} + BB_{2})z_t+AB_{2}z_ts_t + R_tB_{2}z_t + \varepsilon_t}
#'                                                                                                                        \deqn{= \Theta_{0}z_t + \Theta_{1}z_ts_t+\varepsilon_t^{*}}
#'                                                                                                                        where \eqn{\Theta_{0} = B_{1} + B_{2}'B}, \eqn{\Theta_{1} = B_{2}'A} and \eqn{\varepsilon_t^{*} = R_tB_{2} + \varepsilon_t}. Under the null,
#'                                                                                                                        \eqn{\Theta_{0} = B_{1}} and \eqn{\Theta_{1} = 0}, while the previous model is linear, with \eqn{\varepsilon_t^{*} = \varepsilon_t}. It
#'                                                                                                                        follows that the Lagrange multiplier test, under the null, is derived from the score
#'                                                                                                                        \deqn{\frac{\partial \log L(\widetilde{\theta})}{\partial \Theta_{1}} = \sum_{t=1}^{T}z_ts_t(y_t - \widetilde{B}_{1}z_t)'\widetilde{\Omega}^{-1} = S(Y - Z\widetilde{B}_{1})\widetilde{\Omega}^{-1},}
#'                                                                                                                        where
#'                                                                                                                        \deqn{S = z_{1}'s_{1}\\\vdots\\ z_t's_t}
#'                                                                                                                        and where \eqn{\widetilde{B}_{1}} and \eqn{\widetilde{\Omega}} are estimated from the model in \eqn{H_{0}}. If \eqn{P_{Z} = Z(Z'Z)^{-1}Z'} is the
#'                                                                                                                        projection matrix of Z, the LM test is specified as follows
#'                                                                                                                        \deqn{LM = tr\big\{\widetilde{\Omega}^{-1}(Y - Z\widetilde{B}_{1})'S\big[S'(I_t - P_{Z})S\big]^{-1}S'(Y-Z\widetilde{B}_{1})\big\}.}
#'                                                                                                                        Under the null, the test statistics is distributed as a \eqn{\chi^{2}} with \eqn{\widetilde{n}(p\cdot\widetilde{n} + k)} degrees of freedom.
#'
#' @param y \code{data.frame} or \code{matrix} of dependent variables of dimension \code{(Txn)}
#' @param exo (optional) \code{data.frame} or \code{matrix} of exogenous variables of dimension \code{(Txk)}
#' @param st single transition variable for all the equation of dimension \code{(Tx1)}
#' @param st.choice boolean identifying whether the transition variable should be selected from a matrix of \code{R} potential variables of dimension \code{(TxR)}
#' @param alpha Confidence level
#' @return An object of class \code{VLSTARjoint}.
#' @references Luukkonen R., Saikkonen P. and Terasvirta T. (1988), Testing Linearity Against Smooth Transition Autoregressive Models. \emph{Biometrika}, 75: 491-499
#'
#' Terasvirta T. and Yang Y. (2015), Linearity and Misspecification Tests for Vector Smooth Transition Regression Models. \emph{CREATES Research Paper 2014-4}
#' @author Andrea Bucci
#' @aliases print.VLSTARjoint
#' @keywords VLSTAR
#' @export
#' @importFrom vars VAR
#' @importFrom stats lm qchisq pchisq residuals
#' @importFrom matrixcalc matrix.trace
VLSTARjoint <- function(y, exo = NULL, st, st.choice = FALSE, alpha = 0.05){
  y <- as.matrix(y)
  ncoly <- ncol(y)

  x <- exo

##VAR Estimation Y on X
if(!is.null(x)){
  varest <- VAR(y, exogen = x)
  #q <- ncolx-ncoly
}else{
  varest <- VAR(y)
}


  if(alpha > 1 | alpha < 0){
    stop('Please provide a valid value for alpha')
  }

if(ncol(as.matrix(st))>1 & st.choice == FALSE){
  stop('Please provide a valid transition variable or specify the automatic selection')
}

x <- varest$datamat[,-c(1:ncoly)]
ncolx <- ncol(x)
nrowx <- nrow(x)


if(st.choice == TRUE){
  LM3 <- matrix(nrow = ncol(st), ncol = 1)
  pvalue <- matrix(nrow = ncol(st), ncol = 1)
  for(j in 1:ncol(st)){
  st <- as.matrix(st[varest$p:nrow(st),])
  ee <- residuals(varest)
  RSS0 <- t(ee)%*%ee
  ZZ <- matrix(nrow = nrowx, ncol = ncolx*3)
  for (i in 1:nrowx){
    xst1 <-  as.matrix(x[i,]*st[i,j])
    xst2 <- as.matrix(x[i,]*st[i,j]^2)
    xst3 <- as.matrix(x[i,]*st[i,j]^3)
    ZZ[i,] <- cbind(xst1, xst2, xst3)
  }

  ausvar <- lm(ee ~ ZZ)
  ll <- residuals(ausvar)
  RSS1 <- t(ll)%*%ll
  trac1 <- matrix.trace(ginv(RSS0)%*%RSS1)
  LM3[j,] <- nrowx*(ncoly - trac1)
  df <- 3*ncoly + (ncolx)
  conflev <- 1-alpha/2
  chi <- qchisq(conflev, df)
  pvalue[j,] <- pchisq(LM3[j,], df, lower.tail=FALSE)
  }
}else{
  st <- as.matrix(st[varest$p:length(st)])
  ee <- residuals(varest)
  RSS0 <- t(ee)%*%ee
  ZZ <- matrix(nrow = nrowx, ncol = ncolx*3)
  for (i in 1:nrowx){
    xst1 <-  as.matrix(x[i,]*st[i])
    xst2 <- as.matrix(x[i,]*st[i]^2)
    xst3 <- as.matrix(x[i,]*st[i]^3)
    ZZ[i,] <- cbind(xst1, xst2, xst3)
  }
  ausvar <- lm(ee ~ ZZ)
  ll <- residuals(ausvar)
  RSS1 <- t(ll)%*%ll
  trac1 <- matrix.trace(ginv(RSS0)%*%RSS1)
  LM3 <- nrowx*(ncoly - trac1)
  df <- 3*ncoly + (ncolx)
  conflev <- 1-alpha/2
  chi <- qchisq(conflev, df)
  pvalue <- pchisq(LM3, df, lower.tail=FALSE)
}

results <- list(LM3, pvalue, chi, st, st.choice, df)
names(results) <- c('LM', 'pval', 'critical','st', 'st.choice', 'df')
class(results) = 'VLSTARjoint'
#cat('Joint linearity test (Third-order Taylor expansion)\n')
return(results)
}



#' @export
print.VLSTARjoint <- function(x, ...)
{
  if(x$st.choice == TRUE){
    digits = 5
    LM <- max(x$LM)
    cat("\nJoint linearity test (Third-order Taylor expansion)\n")
    cat("Transition variable chosen:", colnames(x$st)[which.max(x$LM)],"\n")
    cat(" LM =", format(LM, digits=digits),"; p-value =", format(x$pval[which.max(x$LM)], digits=digits),"\n")
    cat(" Critical value for alpha =", format(x$critical, digits=digits), "\n")
  }else{
    digits = 3
    cat("\nJoint linearity test (Third-order Taylor expansion)\n")
    cat(" LM =", format(x$LM, digits=digits),"; p-value =", format(x$pval, digits=digits),"\n")
    cat(" Critical value for alpha =", format(x$critical, digits=digits), "\n")
  }
  invisible(x)
}
