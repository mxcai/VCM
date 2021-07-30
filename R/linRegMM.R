# VCM implemented by MM algorithm
linRegMM <- function(X,y,Z=NULL,maxIter=1500,tol=1e-6,se2=NULL,sb2=NULL,verbose=T){
  p <- ncol(X)
  n <- length(y)

  ym <- mean(y)
  y  <- y-ym  # centerize y
  Xm <- colMeans(X)
  X <- scale(X,center = T,scale = F)  # centerize X
  Xsd <- sqrt(colMeans(X^2))
  X <- t(t(X)/Xsd)/sqrt(p)  # scale X
  K <- X %*% t(X)

  if(is.null(Z)){
    Z <- matrix(1,n,1)
  } else{
    Z <- cbind(1,Z)
  }

  #calculate sth in advance
  if(ncol(Z)==1){
    SZy <- 0
  } else {
    ZZ <- t(Z)%*%Z
    SZy <- solve(ZZ,t(Z)%*%y)
  }

  eigenK <- eigen(K)
  eVal <- eigenK$values
  eVec <- eigenK$vectors

  yt <- t(eVec)%*%y
  Zt <- t(eVec)%*%Z

  #initialize
  if(is.null(se2)) {se2 <- drop(var(y))/2}
  if(is.null(sb2)) {sb2 <- drop(var(y))/2}
  lb <- rep(0,maxIter)
  lb[1] <- -Inf

  for(iter in 2:maxIter){

    D <- eVal/se2 + 1/sb2
    d <- 1/(D*se2*sb2) # <=> d <- 1/(sb2*eVal + se2) in the original paper, d = 1/D/se2/sb2
    beta0 <- solve((t(Zt*d))%*%Zt) %*% (t(Zt) %*% (yt*d))
    res <- yt - Zt %*% beta0

    sb2 <- sb2 * sqrt(sum(res^2 * d^2 * eVal) / sum(d * eVal))
    se2 <- se2 * sqrt(sum(res^2 * d^2) / sum(d))
    # safe guard for small variance components
    sb2 <- ifelse(sb2<1e-6,1e-6,sb2)
    se2 <- ifelse(se2<1e-6,1e-6,se2)

    #evaluate lower bound
    lb[iter] <- get_lb_MM(n,d,res)

    if(verbose){
      cat(iter,"-th iteration, lower bound = ",lb[iter]," ,diff=",lb[iter]-lb[iter-1],",sb2=",sb2,",se2=",se2,"\n")
    }
    if(abs(lb[iter]-lb[iter-1])<tol){
      lb <- lb[1:iter]
      break
    }
  }
  mu <- 1/se2 * t(X) %*% (eVec %*% ((yt - Zt%*%beta0) / D))  # after covnergence obtain mu using one E-step
  # recover beta0 and mu
  mu <- mu/Xsd/sqrt(p)
  if(ncol(Z)==1){
    beta0 <- beta0 + ym - colSums(mu*Xm)
  } else{
    beta0 <- beta0 + solve(ZZ,colSums(Z)*(ym-sum(mu*Xm)))  # = beta0 + (Z^TZ)^-1Z^T[(y^bar-X^bar%*%mu),...,(y^bar-X^bar%*%mu)]
  }


  invSigy <- eVec%*%(1/(eVal*sb2+se2)*t(eVec))    # solve(sb2*K+se2*diag(n))
  invSigyK <- invSigy%*%K
  FIM <- matrix(0,2,2)
  FIM[1,1] <- sum(invSigyK^2) / 2
  FIM[2,2] <- sum(invSigy^2) / 2
  FIM[1,2] <- FIM[2,1] <- sum(invSigyK*invSigy) / 2
  covSig <- solve(FIM)  #inverse of FIM

  bayesRegMM <- list(beta0=beta0,sb2=sb2,se2=se2,mu=mu,K=K,iter=iter,covSig=covSig,lb=lb)

  class(bayesRegMM) <- c("VCM","MM")
  bayesRegMM
}


####################################################################################################
get_lb_MM <- function(n,d,res){
  llh <- sum(log(d))/2 - sum(res^2*d)/2 - n/2*log(2*pi)
}
