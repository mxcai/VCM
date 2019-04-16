predict.mvVCM <- function(object,Xnew,Znew=NULL){
  nY <- length(object$sb2)
  n <- nrow(Xnew)
  q <- nrow(object$beta0)

  if(n!=nrow(Znew)) stop("Inputs X and Z do not match.")

  if(is.null(Znew)) {
    Znew <- matrix(1,n,1)
  } else {
    Znew <- cbind(1,Znew)
  }

  if(q!=ncol(Znew)) stop("Inputs arguments object and Z are not compatible.")

  Y_hat <- matrix(0,n,nY)
  for(k in 1:nY){
    Y_hat[,k] <- Znew%*%object$beta0[,k] + Xnew%*%object$mu[,k]
  }

  return(Y_hat)
}


predict.VCM <- function(object,Xnew,Znew=NULL){
  n <- nrow(Xnew)
  q <- nrow(object$beta0)

  if(n!=nrow(Znew)) stop("Inputs X and Z do not match.")

  if(is.null(Znew)) {
    Znew <- matrix(1,n,1)
  } else {
    Znew <- cbind(1,Znew)
  }

  if(q!=ncol(Znew)) stop("Inputs arguments object and Z are not compatible.")

  y_hat <- Znew%*%object$beta0 + Xnew%*%object$mu

  return(y_hat)
}
