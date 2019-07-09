multiCUMSUM <- function(data, realized){
  nrowy <- nrow(data)
  ncoly <- ncol(data)
  dmat <- as.matrix(data)
  Sk <- matrix(nrow = nrowy, ncol = ncoly)

for(k in 1 : nrowy){
  Sk0 <- vech(dmat[k,]%*%t(dmat[k,])) + Sk0
  Sk[k,] <- (1/nrowy^0.5)*(Sk0 - (k/nrowy)*Sk1)
}

}

