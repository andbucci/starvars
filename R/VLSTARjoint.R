VLSTARjoint <- function(y, exo = NULL, st, alpha = 0.05){
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

x <- varest$datamat[,-c(1:ncoly)]
ncolx <- ncol(x)
nrowx <- nrow(x)
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
df <- 3*ncoly*(ncolx)
conflev <- 1-alpha/2
chi <- qchisq(conflev, df)
pvalue <- pchisq(LM3, df, lower.tail=FALSE)
results <- list(LM3, pvalue, chi)
names(results) <- c('LM', 'pval', 'critical')
class(results) = 'VLSTARjoint'
#cat('Joint linearity test (Third-order Taylor expansion)\n')
return(results)
}



#' @S3method print.VLSTARjoint
print.VLSTARjoint <- function(x, ...)
{
  digits = 3
  cat("\nJoint linearity test (Third-order Taylor expansion)\n")
  cat(" LM =", format(x$LM, digits=digits),"; p-value =", format(x$pval, digits=digits),"\n")
  cat(" Critical value for alpha =", format(x$critical, digits=digits), "\n")
  invisible(x)
}
