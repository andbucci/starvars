multiCUMSUM <- function(data, conf.level = 0.95){

ncoly <- ncol(data)
nrowy <- nrow(data)

vechy <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
for (j in 1:nrowy){
  vechy[j,] <- vech(t(data[j,])%*%data[j,])
}

sumvechy <- colSums(vechy)

buildtau <- function(data){
  V <- as.matrix(data)
  T1 <- nrowy
  nupper <- 0.5*ncoly*(ncoly+1)
  tau <- matrix(0L, nrow = T1, ncol = nupper)
  for(t in 1:T1){
    tau[t,] <- vech(V[t,]%*%t(V[t,]))
  }
  zeros <- as.matrix(t(rep(0, nupper)))
  tau <- rbind(zeros, tau)
  tau <- zoo(tau, order.by = index(data))
  return(tau)
}


makeDT <- function(tau){
  nupper <- ncol(tau)
  DT <- matrix(0L, ncol = nupper, nrow = nupper)
  for(i in 1:nupper){
    DT[i, i] <- starvars::lrvarbart(tau[,i])$lrv
  }
  return(DT)
}

tau <- buildtau(data)
D <- makeDT(tau)

Sk <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
tmp1 <- matrix(0, ncol = 1, nrow = (ncoly*(ncoly+1)/2))
Lambda_hat <- NULL
for (i in 1:nrowy){
  tmp <- vech(t(data[i,])%*%data[i,])
    tmp1 <- tmp1 + tmp
    Sk[i,] <- t(1/sqrt(nrow(data)) * (tmp1 - (i/nrowy)*sumvechy))
    Lambda_hat[i] <- t(Sk[i,])%*%ginv(D)%*%Sk[i,]

}

Lambda <- max(Lambda_hat)
Omega <- mean(Lambda_hat)

breakpoint <- as.Date(index(data)[which.max(Lambda_hat)], origin="1970-01-01")


critLambda <- rbind(c(2.64,2.53,2.46,2.27,2.16,2.06,1.96,1.28),
                    c(3.17,3.02,2.92,2.69,2.55,2.44,2.33,1.64),
                    c(4.28,4.04,3.89,3.53,3.33,3.18,3.04,2.33))
critOmega <- rbind(c(1.33,1.33,1.32,1.31,1.31,1.30,1.29,1.28),
                   c(1.84,1.81,1.79,1.74,1.71,1.69,1.68,1.64),
                   c(2.90,2.80,2.74,2.59,2.51,2.46,2.41,2.33))

d1 <- ncol(tau)

if(conf.level == 0.9){
  r = 1
}else if(conf.level == 0.95){
  r = 2
}else if(conf.level == 0.99){
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

M = cbind(Lambda, Omega)
breaklocation <- data.frame(which.max(Lambda_hat), breakpoint)
colnames(breaklocation) <- c('index', 'Date')
multiCS <- list(M,breaklocation, r,  c1, critLambda, critOmega)
names(multiCS) <- c('M', 'break_location', 'r', 'c', 'critLambda', 'critOmega')
class(multiCS) <- 'multiCUMSUM1'
return(multiCS)
}


print.multiCUMSUM <- function(x, ...) {
  cat("===============================================\n")
  cat("Break detection in the covariance structure:\n")
  cat('Lambda (d) test statistics: ',x$M[1], ' [', x$critLambda[x$r, x$c], ']\n', sep = '')
  cat('Omega (d) test statistics: ', x$M[2], ' [', x$critOmega[x$r, x$c], ']\n', sep = '')
  if(x$M[1] > x$critLambda[x$r, x$c] | x$M[2] > x$critOmega[x$r, x$c]){
    cat('Break located at Date: ', as.character(x$break_location$Date), sep = '')
  }else{
    cat('No common breaks were found.')
  }

}
