#'@S3method predict VLSTAR

predict.VLSTAR <- function(object, newdata, alpha = 0.05, ...){
  k <- ncol(object$Data[[2]])
  newdata <- as.matrix(c(object$Data[[1]][(nrow(object$Data[[1]])-object$p):(nrow(object$Data[[1]])-object$p),], as.matrix(newdata)))
  if(!missing(newdata)) {
    if(!inherits(newdata, c("data.frame", "matrix","zoo", "ts"))) stop("Arg 'newdata' should be of class data.frame, matrix, zoo or ts")
    if(!dim(newdata)== (k-1))stop("newdata should have the same dimension of the design matrix (p*n+k)")}
  ## extract parameters, coefs
  if(object$constant == T){
    newdata <- rbind(1, newdata)
  }
  newdata <- as.matrix(newdata)

  BB <- object$B
  In <- diag(ncol(object$Data[[1]])
  Gtilde <- t(cbind(In, object$Gtilde))
  pred <- as.data.frame((t(Gtilde)%*%t(BB)%*%newdata))
  rownames(pred) <- colnames(object$Bhat)
  colnames(pred) <- 'Prediction'
  pred$LB <- NA
  pred$UB <- NA
  for(i in 1:ncol(object$Data[[1]])){
   std_err <- sqrt(object$Omega[i,i])*t(newdata)%*%MASS::ginv(t(object$Data[[2]])%*%object$Data[[2]])%*%newdata
     pred$LB[i] <- pred$Prediction[i]-stats::qt((1-alpha), df = nrow(object$Data[[1]])*ncol(object$Data[[1]])-length(BB))*std_err
     pred$UB[i] <- pred$Prediction[i]+stats::qt((1-alpha), df = nrow(object$Data[[1]])*ncol(object$Data[[1]])-length(BB))*std_err
  }

  return(pred)
}
