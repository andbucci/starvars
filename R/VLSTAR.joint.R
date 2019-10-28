VLSTAR.joint <- function(y1, x1, st, alpha = 0.05){
  require(matrixcalc)
  y <- as.matrix(y1)
  x <- as.matrix(x1)
  nrowx <- nrow(x)
  ncolx <- ncol(x)
  ncoly <- ncol(y)
  q <- ncolx-ncoly

  ##VAR Estimation Y on X
varest <- lm(y~x)
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
df <- 3*ncoly*(q+ncoly)
conflev <- 1-alpha/2
chi <- qchisq(conflev, df)
pvalue <- pchisq(LM3, df, lower.tail=FALSE)
results <- list(LM3, pvalue, chi)
names(results) <- c('LM', 'p-value', 'Critical value for alpha')
cat('Joint linearity test (Third-order Taylor expansion)\n')
return(results)
}
