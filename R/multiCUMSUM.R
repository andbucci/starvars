multiCUMSUM <- function(data, conf.level = 0.95, max.breaks = 7){

  if(conf.level > 1 | conf.level < 0){
    stop('Please provide a valid value for alpha')
  }
nrowy <- nrow(data)
ncoly <- ncol(data)

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
  tau <- as.xts(tau, order.by = index(data))
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

critLambda <- rbind(c(2.64,2.53,2.46,2.27,2.16,2.06,1.96,1.28),
                    c(3.17,3.02,2.92,2.69,2.55,2.44,2.33,1.64),
                    c(4.28,4.04,3.89,3.53,3.33,3.18,3.04,2.33))
critOmega <- rbind(c(1.33,1.33,1.32,1.31,1.31,1.30,1.29,1.28),
                   c(1.84,1.81,1.79,1.74,1.71,1.69,1.68,1.64),
                   c(2.90,2.80,2.74,2.59,2.51,2.46,2.41,2.33))


Lambda <- list()
Omega <- list()
breakpoint <- list()
vechy <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
for (j in 1:nrowy){
  vechy[j,] <- vech(t(data[j,])%*%data[j,])
}
sumvechy <- colSums(vechy)
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
Lambda[[1]] <- max(Lambda_hat)
Omega[[1]] <- mean(Lambda_hat)
breakpoint[[1]] <- as.Date(index(data)[which.max(Lambda_hat)], origin="1970-01-01")
d1 <- ncol(tau)

ifelse(nrow(data[1:which.max(Lambda_hat),])==1, data1 <- data[1:which.max(Lambda_hat),], data1 <- data[1:(which.max(Lambda_hat)-1),])
ifelse(nrow(data[which.max(Lambda_hat):nrow(data),])==1, data2 <- data[which.max(Lambda_hat):nrow(data),], data2 <- data[(which.max(Lambda_hat)+1):nrow(data),])
data_1 <- list(data, data1, data2)

for(k in 2:3){
  nrowy <- nrow(data_1[[k]])
  if(nrowy == 1){
    Lambda[[k]] <- NA
    Omega[[k]] <- NA
    breakpoint[[k]] <- NA
  }else{
    vechy <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
    for (j in 1:nrowy){
      vechy[j,] <- vech(t(data_1[[k]][j,])%*%data_1[[k]][j,])
    }
    sumvechy <- colSums(vechy)
    tau <- buildtau(data_1[[k]])
    D <- makeDT(tau)
    Sk <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
    tmp1 <- matrix(0, ncol = 1, nrow = (ncoly*(ncoly+1)/2))
    Lambda_hat <- NULL
    for (i in 1:nrowy){
      tmp <- vech(t(data_1[[k]][i,])%*%data_1[[k]][i,])
      tmp1 <- tmp1 + tmp
      Sk[i,] <- t(1/sqrt(nrow(data_1[[k]])) * (tmp1 - (i/nrowy)*sumvechy))
      Lambda_hat[i] <- t(Sk[i,])%*%ginv(D)%*%Sk[i,]

    }
    Lambda[[k]] <- max(Lambda_hat)
    Omega[[k]] <- mean(Lambda_hat)
    breakpoint[[k]] <- as.Date(index(data_1[[k]])[which.max(Lambda_hat)], origin="1970-01-01")
  }
}

ifelse(nrow(data1[1:which(index(data1)==breakpoint[[2]]),])==1, data3 <- data1[1:which(index(data1)==breakpoint[[2]]),], data3 <- data1[1:(which(index(data1)==breakpoint[[2]])-1),])
ifelse(nrow(data1[which(index(data1)==breakpoint[[2]]):nrow(data1),])==1, data4 <- data1[which(index(data1)==breakpoint[[2]]):nrow(data1),], data4 <- data1[(which(index(data1)==breakpoint[[2]])+1):nrow(data1),])
ifelse(nrow(data2[1:which(index(data2)==breakpoint[[3]]),])==1, data5 <- data2[1:which(index(data2)==breakpoint[[3]]),], data5 <- data2[1:(which(index(data2)==breakpoint[[3]])-1),])
ifelse(nrow(data2[which(index(data2)==breakpoint[[3]]):nrow(data2),])==1, data6 <- data2[which(index(data2)==breakpoint[[3]]):nrow(data2),], data6 <- data2[(which(index(data2)==breakpoint[[3]])+1):nrow(data2),])


data_1 <- list(data, data1, data2, data3, data4, data5, data6)

for(k in 4:7){
  nrowy <- nrow(data_1[[k]])
  if(nrowy == 1){
    Lambda[[k]] <- NA
    Omega[[k]] <- NA
    breakpoint[[k]] <- NA
  }else{
    vechy <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
    for (j in 1:nrowy){
      vechy[j,] <- vech(t(data_1[[k]][j,])%*%data_1[[k]][j,])
    }
    sumvechy <- colSums(vechy)
    tau <- buildtau(data_1[[k]])
    D <- makeDT(tau)
    Sk <- matrix(ncol = (ncoly*(ncoly+1)/2), nrow = nrowy)
    tmp1 <- matrix(0, ncol = 1, nrow = (ncoly*(ncoly+1)/2))
    Lambda_hat <- NULL
    for (i in 1:nrowy){
      tmp <- vech(t(data_1[[k]][i,])%*%data_1[[k]][i,])
      tmp1 <- tmp1 + tmp
      Sk[i,] <- t(1/sqrt(nrow(data_1[[k]])) * (tmp1 - (i/nrowy)*sumvechy))
      Lambda_hat[i] <- t(Sk[i,])%*%ginv(D)%*%Sk[i,]

    }
    Lambda[[k]] <- max(Lambda_hat)
    Omega[[k]] <- mean(Lambda_hat)
    breakpoint[[k]] <- as.Date(index(data_1[[k]])[which.max(Lambda_hat)], origin="1970-01-01")
  }
}
Lambda <- unlist(Lambda)
breakpoint <- unlist(breakpoint)
breaks <- list()
breaks[[1]] <- as.Date(breakpoint[1])
breaks[[2]] <- c(as.Date(breakpoint[1]), as.Date(breakpoint[2:3][which.max(Lambda[2:3])]))
breaks[[3]] <- c(as.Date(breakpoint[1]), as.Date(breakpoint[2:3][which.max(Lambda[2:3])]), as.Date(breakpoint[2:3][which.min(Lambda[2:3])]))
breaks[[4]] <- c(as.Date(breakpoint[1]), as.Date(breakpoint[2:3][which.max(Lambda[2:3])]), as.Date(breakpoint[2:3][which.min(Lambda[2:3])]),
                 as.Date(breakpoint[4:7][which.max(Lambda[4:7])]))
breaks[[5]] <- c(as.Date(breakpoint[1]), as.Date(breakpoint[2:3][which.max(Lambda[2:3])]), as.Date(breakpoint[2:3][which.min(Lambda[2:3])]),
                 as.Date(breakpoint[4:7][which.max(Lambda[4:7])]), as.Date(breakpoint[which(Lambda == sort(Lambda[4:7])[3])]))
breaks[[6]] <- c(as.Date(breakpoint[1]), as.Date(breakpoint[2:3][which.max(Lambda[2:3])]), as.Date(breakpoint[2:3][which.min(Lambda[2:3])]),
                 as.Date(breakpoint[4:7][which.max(Lambda[4:7])]), as.Date(breakpoint[which(Lambda == sort(Lambda[4:7])[3])]),
                 as.Date(breakpoint[which(Lambda == sort(Lambda[4:7])[2])]))
breaks[[7]] <- c(as.Date(breakpoint[1]), as.Date(breakpoint[2:3][which.max(Lambda[2:3])]), as.Date(breakpoint[2:3][which.min(Lambda[2:3])]),
                 as.Date(breakpoint[4:7][which.max(Lambda[4:7])]), as.Date(breakpoint[which(Lambda == sort(Lambda[4:7])[3])]),
                 as.Date(breakpoint[which(Lambda == sort(Lambda[4:7])[2])]), as.Date(breakpoint[which(Lambda == sort(Lambda[4:7])[1])]))


Lambda1 <- matrix(ncol = 1, nrow = 7)
Lambda1[1] <- Lambda[1]
Lambda1[2] <- max(Lambda[2:3])
Lambda1[3] <- min(Lambda[2:3])
Lambda1[4] <- sort(Lambda[4:7])[4]
Lambda1[5] <- sort(Lambda[4:7])[3]
Lambda1[6] <- sort(Lambda[4:7])[2]
Lambda1[7] <- sort(Lambda[4:7])[1]

Omega <- unlist(Omega)
Omega1 <- matrix(ncol = 1, nrow = 7)
Omega1[1] <- Omega[1]
Omega1[2] <- max(Omega[2:3])
Omega1[3] <- min(Omega[2:3])
Omega1[4] <- sort(Omega[4:7])[4]
Omega1[5] <- sort(Omega[4:7])[3]
Omega1[6] <- sort(Omega[4:7])[2]
Omega1[7] <- sort(Omega[4:7])[1]

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
M <- M[1:max.breaks,]
breaks <- breaks[1:max.breaks]
multiCS <- list(M,breaks, r,  c1, critLambda, critOmega)
names(multiCS) <- c('M', 'break_location', 'r', 'c', 'critLambda', 'critOmega')
class(multiCS) <- 'multiCUMSUM'
return(multiCS)
}


print.multiCUMSUM <- function(x, ...) {
  nrowlam <- nrow(x$M)
  breakpoint <- matrix(ncol = nrowlam, nrow = nrowlam)
  for(i in 1:nrowlam){
    breakpoint[i,] <- c(as.character(x$break_location[[i]]), rep('', (nrowlam-i)))
  }
  MM <- as.data.frame(round(x$M,2))
  MM <- cbind(MM, breakpoint)

  row.names(MM) <- paste('Break ', 1:nrowlam, sep ='')
  colnames(MM)[3:(nrowlam+2)] <- paste('Break Date ', 1:nrowlam, sep ='')
  cat(rep("============", (nrowlam+2)), sep ='', "\n")
  cat("Break detection in the covariance structure:\n")
  print(MM)
  cat(rep("============", (nrowlam+2)), sep ='', "\n")
  cat("Critical values are" , x$critLambda[x$r, x$c], "for Lambda and", x$critOmega[x$r, x$c], "for Omega.\n")
  if(any(x$M[,1] > x$critLambda[x$r, x$c]) | any(x$M[,2] > x$critOmega[x$r, x$c])){
    cat(which.max(x$M[,1]), 'Break(s) identified with Lambda', '\n')
    cat(which.max(x$M[,2]), 'Break(s) identified with Omega', '\n')
  }else{
    cat('No common breaks were found.')
  }

}
