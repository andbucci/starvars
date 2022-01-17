#' Summary method for objects of class VLSTAR
#'
#' \sQuote{\code{summary}} methods for class \sQuote{\code{VLSTAR}}.
#' @aliases summary print.summary.VLSTAR print.summary
#' @param object An object of class \sQuote{\code{VLSTAR}} obtained through \command{VLSTAR()}.
#' @param x A summary object of class \sQuote{\code{VLSTAR}} obtained through \command{summary()}.
#' @param \dots further arguments to be passed to and from other methods
#' @references Terasvirta T. and Yang Y. (2014), Specification, Estimation and Evaluation of Vector Smooth Transition Autoregressive Models with Applications. \emph{CREATES Research Paper 2014-8}
#' @author Andrea Bucci
#' @return An object of class \code{summary.VLSTAR} containing a list of summary information from VLSTAR estimates. When \code{print} is applied to this object, summary information are printed
#' @keywords VLSTAR
#' @seealso \code{\link{VLSTAR}}
#' @export
#'
summary.VLSTAR<-function(object,...){
  x<-object
  k<-ncol(x$Data[[2]])
  t<-nrow(x$Data[[1]])
  p<-x$p
  x$T <- nrow(x$yoriginal)
  Z<-t(as.matrix(tail.matrix(x$Data[[1]])))
  x$npar <- 2*ncol(x$Data[[1]])*x$m + x$m*ncol(x$Data[[2]])*ncol(x$Data[[1]])

  ## export results
  x$coefficients<-as.list(as.data.frame(x$Bhat))
  x$StDev<-as.list(as.data.frame(x$StDev))
  x$Pvalues<-as.list(as.data.frame(x$pval))
  x$Tvalues<-as.list(as.data.frame(x$ttest))
  ab<-list()
  symp<-list()
  stars<-list()
  for(i in 1:length(x$Pvalues)){
    symp[[i]] <- symnum(x$Pvalues[[i]], corr=FALSE,cutpoints = c(0,  .001,.01,.05, .1, 1), symbols = c("***","**","*","."," "))
    stars[[i]]<-matrix(symp[[i]], nrow=length(x$Pvalues[[i]]))
    ab[[i]]<-matrix(paste(x$coefficients[[i]],"(", x$StDev[[i]],")",stars[[i]], sep=""), nrow=length(x$StDev[[i]]))
    dimnames(ab[[i]])<-dimnames(x$coefficients[[1]])
  }
  attributes(ab)<-attributes(x$coefficients)
  x$stars<-stars
  x$starslegend<-symp[[1]]
  x$bigcoefficients<-ab
  x$aic<-x$AIC
  x$bic<-x$BIC
  x$t <- t
  x$k <- k
  class(x) <- 'summary.VLSTAR'
  return(x)
  NextMethod('print.summary')
}

#' @export
#' @describeIn summary.VLSTAR Print of the summary
print.summary.VLSTAR<-function(x,...){
  digits = 3
  coeftoprint<-list()
  myformat<-function(x,digits, toLatex=FALSE){
    r<-x
    littlex<-abs(x)<10^-(digits)
    r[!littlex]<-formatC(x[!littlex],digits=digits, format="f")
    r[littlex]<-format(x[littlex],digits=min(digits,2), scientific=TRUE)
    if(toLatex)
      r<-gsub("(e.*)","slashtext{\\1}",r)
    if (is(x, "numeric"))
      return(noquote(r))
    if (is(x, "matrix"))
      return(matrix(noquote(r), ncol=ncol(x), nrow=nrow(x)))
  }
  for(i in 1:length(x$bigcoefficients)){
    a<-myformat(x$coefficients[[i]], digits)
    #b<-myformat(x$StDev[[i]], digits)
    #c<-myformat(x$Pvalues[[i]], digits)
    aic1 <- round(x$AIC,2)
    bic1 <- round(x$BIC,2)
    ll1 <- round(x$LL, 2)
    stars1<-x$stars[[i]]
    coeftop <- cbind(a, stars1)
    coeftop1 <- rep(NA, )
    for(j in 1:nrow(coeftop)){
      coeftop1[j] <- paste(coeftop[j,1], coeftop[j,2], sep = '')
    }
    names(coeftop1) = rownames(x$Bhat)
    #colnames(coeftop) <- c('Estimate','stars')
    if(x$singlecgamma == TRUE){
      coeftoprint[[i]] <- list(coefficients = coeftop1, AIC = aic1[i], BIC = bic1[i], LL = ll1[i])
    } else{
      coeftoprint[[i]] <- list(coefficients = coeftop1, gamma =  round(x$Gammac[seq(i, nrow(x$Gammac), ncol(x$Data[[1]])),1],4),
                               c = round(x$Gammac[seq(i, nrow(x$Gammac), ncol(x$Data[[1]])),2],4),
                               AIC = aic1[i], BIC = bic1[i], LL = ll1[i])
    }

  }
  names(coeftoprint) <- colnames(x$Bhat)
  cat("Model VLSTAR with ", x$m, " regimes\n", sep ='')
  cat("Full sample size:",x$T, "\n")
  cat("Number of estimated parameters:", x$npar, "\tMultivariate log-likelihood:", x$MultiLL,"\n")

  if(x$singlecgamma == TRUE){
    cat("\nUnique gamma:", round(x$Gammac[seq(1, nrow(x$Gammac), ncol(x$Data[[1]])),1],4),"\tUnique c:", round(x$Gammac[seq(1, nrow(x$Gammac), ncol(x$Data[[1]])),2],4), "\n")
  }
  cat("==================================================\n")
  #cat('\nEquation:\n')
  #print(noquote(coeftoprint))
  for(i in 1:length(x$bigcoefficients)){
    cat("\nEquation", paste("y", i, sep =''), "\n")

    for(k in 1:x$m){
      cat("\nCoefficients regime", k, "\n")
      coeftmp = coeftoprint[[i]]$coefficients[grepl(paste("m_", k, sep = ''), names(coeftoprint[[i]]$coefficients))]
      names(coeftmp) = colnames(x$Data[[2]])
      print(noquote(coeftmp))
    }
    if(x$singlecgamma == FALSE){
    if(x$m == 2){
      cat("\nGamma:", coeftoprint[[i]]$gamma, "\tc:", coeftoprint[[i]]$c, "\n")
    }else{
      cat("\nGamma:", coeftoprint[[i]]$gamma, "\n")
      cat("\nc:", coeftoprint[[i]]$c, "\n")
    }
      }


    cat("AIC:", coeftoprint[[i]]$AIC, "\tBIC:", coeftoprint[[i]]$BIC, "\tLL:", coeftoprint[[i]]$LL, "\n")

    #cat("==================================================\n")
  }
  cat("==================================================\n")
  cat("\n Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
}
