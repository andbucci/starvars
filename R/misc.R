myformat<-function(x,digits, toLatex=FALSE){
  r<-x
  littlex<-abs(x)<10^-(digits)
  r[!littlex]<-formatC(x[!littlex],digits=digits, format="f")
  r[littlex]<-format(x[littlex],digits=min(digits,2), scientific=TRUE)
  if(toLatex)
    r<-gsub("(e.*)","slashtext{\\1}",r)
  if(class(x)=="numeric")
    return(noquote(r))
  if(class(x)=="matrix")
    return(matrix(noquote(r), ncol=ncol(x), nrow=nrow(x)))
}

###Assign to class list is class matrix
asListIfMat<-function(x){
  if(class(x)=="matrix")
    return(list(x))
  else
    return(x)
}


myInsertCol<-function (m, c, v = NA) {
  #m: matrix
  #c: columns numbers
  #v: value to add (scalar only)
  nr <- nrow(m)
  nc <- ncol(m)
  #first
  m2 <- if (1 %in% c) cbind(matrix(v, nrow = nr), m) else m
  #inter
  for(i in c[!c%in%c(1, length(c)+nc)])
    m2 <- cbind(m2[, 1:(i - 1), drop=FALSE], matrix(v, nrow = nr), m2[,i:ncol(m2), drop=FALSE])
  #last
  if (eval(ncol(m2) + 1) %in% c) 
    m2 <- cbind(m2, matrix(v, nrow = nr))
  
  return(m2)
}
