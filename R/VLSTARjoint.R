VLSTARjoint <- function(y1, x1, st, alpha = 0.05){
  y <- as.matrix(y1)
  x <- as.matrix(x1)
  nrowx <- nrow(x)
  ncolx <- ncol(x)
  ncoly <- ncol(y)
  q <- ncolx-ncoly

  ##VAR Estimation Y on X
varest <- stats::lm(y~x)
ee <- stats::residuals(varest)
RSS0 <- t(ee)%*%ee
ZZ <- matrix(nrow = nrowx, ncol = ncolx*3)
for (i in 1:nrowx){
  xst1 <-  as.matrix(x[i,]*st[i])
  xst2 <- as.matrix(x[i,]*st[i]^2)
  xst3 <- as.matrix(x[i,]*st[i]^3)
 ZZ[i,] <- cbind(xst1, xst2, xst3)
}
ausvar <- stats::lm(ee ~ ZZ)
ll <- stats::residuals(ausvar)
RSS1 <- t(ll)%*%ll
trac1 <- matrixcalc::matrix.trace(MASS::ginv(RSS0)%*%RSS1)
LM3 <- nrowx*(ncoly - trac1)
df <- 3*ncoly*(q+ncoly)
conflev <- 1-alpha/2
chi <- stats::qchisq(conflev, df)
pvalue <- stats::pchisq(LM3, df, lower.tail=FALSE)
results <- list(LM3, pvalue, chi)
names(results) <- c('LM', 'pval', 'critical')
class(results) = 'VLSTARjoint'
#cat('Joint linearity test (Third-order Taylor expansion)\n')
return(results)
}



#' @S3method print.VLSTARjoint
print.VLSTARjoint <- function(object, digits=max(3, getOption("digits") - 2), ...)
{
  x <- object
  cat("\nJoint linearity test (Third-order Taylor expansion)\n")
  cat(" LM =", format(x$LM, digits=digits),"; p-value =", format(x$pval, digits=digits),"\n")
  cat(" Critical value for alpha =", format(x$critical, digits=digits), "\n")
  invisible(x)
}
