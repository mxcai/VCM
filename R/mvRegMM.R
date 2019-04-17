# VCM implemented by MM algorithm
mvRegMM <- function(X,Y,Z=NULL,maxIter=1500,tol=1e-6,se2_init=NULL,sb2_init=NULL,verbose=T){
  Y <- as.matrix(Y)
  nY <- ncol(Y)
  p <- ncol(X)
  n <- nrow(y)

  Ym <- colMeans(Y)
  Y <- scale(Y,center = T,scale = F)  # centerize y
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
    SZY <- rep(0,nY)
  } else {
    ZZ <- t(Z)%*%Z
    SZY <- solve(ZZ,t(Z)%*%Y)
  }

  eigenK <- eigen(K)
  eVal <- eigenK$values
  eVec <- eigenK$vectors

  Yt <- t(eVec)%*%Y
  Zt <- t(eVec)%*%Z

  beta0_all <- matrix(0,ncol(Z),nY)
  sb2_all <- rep(0,nY)
  se2_all <- rep(0,nY)
  mu_all <- matrix(0,p,nY)
  iter_all  <- rep(0,nY)
  covSig_all <- vector("list",nY)
  lb_all <- vector("list",nY)

  for(k in 1:nY){
    y <- Y[,k]
    ym <- Ym[k]
    yt <- Yt[,k]
    if(ncol(Z)==1){
      SZy <- SZY[k]
    } else {
      SZy <- SZY[,k]
    }

    #initialize
    if(is.null(se2_init)) {
      se2 <- drop(var(y))/2
    } else {
      se2 <- se2_init[k]
    }
    if(is.null(sb2_init)) {
      sb2 <- drop(var(y))/2
    } else {
      sb2 <- sb2_init[k]
    }
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
      beta0 <- beta0 + solve(ZZ,colSums(Z)*(ym-sum(mu*Xm)))  # = beta0 - (Z^TZ)^-1Z^T[(y^bar-X^bar%*%mu),...,(y^bar-X^bar%*%mu)]
    }


    invSigy <- eVec%*%(1/(eVal*sb2+se2)*t(eVec))    # solve(sb2*K+se2*diag(n))
    invSigyK <- invSigy%*%K
    FIM <- matrix(0,2,2)
    FIM[1,1] <- sum(invSigyK^2) / 2
    FIM[2,2] <- sum(invSigy^2) / 2
    FIM[1,2] <- FIM[2,1] <- sum(invSigyK*invSigy) / 2
    covSig <- solve(FIM)  #inverse of FIM

    beta0_all[,k] <- beta0
    sb2_all[k] <- sb2
    se2_all[k] <- se2
    mu_all[,k] <- mu
    iter_all[k]  <- iter
    covSig_all[[k]] <- covSig
    lb_all[[k]] <- lb
  }



  bayesRegMM <- list(beta0=beta0_all,sb2=sb2_all,se2=se2_all,mu=mu_all,K=K,iter=iter_all,covSig=covSig_all,lb=lb_all)

  class(bayesRegMM) <- c("mvVCM","MM")
  bayesRegMM
}


####################################################################################################
get_lb_MM <- function(n,d,res){
  llh <- sum(log(d))/2 - sum(res^2*d)/2 - n/2*log(2*pi)
}
