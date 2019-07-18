multiCUMSUM <- function(data, alpha){
  require(ks)
  require(MASS)
  require(sandwich)
  require(fGarch)
  require(R.utils)
nrowy <- nrow(data)+1
ncoly <- ncol(data)
h1 <- matrix(ncol = ncoly, nrow = (nrowy-1))
means <- colMeans(data)
quiet <- function(x){
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
for(i in 1:ncoly){
  gspec <- quiet(garchFit(~garch(1,1), data = data[,i]))
  h1[,i] <- gspec@h.t
  for(j in 1:(nrowy-1)){
  data[j, i] <- (data[j,i] - means[i])/sqrt(h1[j,i])
  }
}

buildtau <- function(data){
    V <- as.matrix(data)
    T1 <- nrowy - 1
    nupper <- 0.5*ncoly*(ncoly+1)
    tau <- matrix(0L, nrow = T1, ncol = nupper)
for(t in 1:T1){
  tau[t,] <- vech(V[t,]%*%t(V[t,]))
}
zeros <- as.matrix(t(rep(0, nupper)))
tau <- rbind(zeros, tau)
return(tau)
}

makeDT <- function(tau){
    nupper <- ncol(tau)
    DT <- matrix(0L, ncol = nupper, nrow = nupper)
    for(i in 1:nupper){
      DT[i, i] <- lrvarbart(tau[,i])$lrv
      }
  return(DT)
}

tau <- buildtau(data)
D <- makeDT(tau)

T1 <- nrow(tau)
d1 <- ncol(tau)
tauT <- colSums(tau)/T1
critLambda <- rbind(c(2.64,2.53,2.46,2.27,2.16,2.06,1.96,1.28),
  c(3.17,3.02,2.92,2.69,2.55,2.44,2.33,1.64),
  c(4.28,4.04,3.89,3.53,3.33,3.18,3.04,2.33))
critOmega <- rbind(c(1.33,1.33,1.32,1.31,1.31,1.30,1.29,1.28),
    c(1.84,1.81,1.79,1.74,1.71,1.69,1.68,1.64),
    c(2.90,2.80,2.74,2.59,2.51,2.46,2.41,2.33))
  if(alpha == 0.9){
    r = 1
  }else if(alpha == 0.95){
    r = 2
  }else if(alpha == 0.99){
    r = 3
  }
  if(d1 <= 10){
    c1 = 1
  }else if(d1<=15 & d1>10){
    c1 = 2
  }else if(d1<= 20 & d1>15){
    c1 = 3
  }else if(d1<=50 & d1>20){
    c1 = 4
  }else if(d1<=100 & d1>50){
    c1 = 5
  }else if(d1<= 200 & d1>100){
    c1 = 6
  }else if(d1<= 500 & d1>200){
    c1 = 7
  }else if(d1 > 500){
    c1 = 8
  }
  M1 = 0
  M2 = 0
  for(t in 1:T1){
    tmp = t(tau[t,]-t*tauT)%*%ginv(D)%*%(tau[t,]-t*tauT)
    if(tmp > M1){
      M1 = tmp
    }
M2 = tmp+M2
  }
M1 = M1/T1
M2 = M2/T1^2

M1 = (M1-d1/4)/sqrt(d1/8)
M2 = (M2-d1/6)/sqrt(d1/45)
M = cbind(M1, M2)

cat("===============================================\n")
cat("Break detection in the covariance structure:\n")
cat('Lambda (d) test statistics: ', M1, ' [', critLambda[r, c1], ']\n', sep = '')
cat('Omega (d) test statistics: ', M2, ' [', critOmega[r, c1], ']\n', sep = '')



return(M)
}
