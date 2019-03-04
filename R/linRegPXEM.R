# ridge regression implemented by PXEM
linRegPXEM <- function(X,y,Z=NULL,maxIter=1500,tol=1e-6,se2=NULL,sb2=NULL,verbose=T){
  p <- ncol(X)
  n <- length(y)

  ym <- mean(y)
  Xm <- colMeans(X)
  X <- scale(X)/sqrt(p)
  K <- X %*% t(X)

  if(is.null(Z)){
    Z <- matrix(1,n,1)
  } else{
    Z <- cbind(1,Z)
  }

  #calculate sth in advance
  if(ncol(Z)==1){
    SZy <- ym
    SZX <- Xm
    # X <- scale(X,scale = F,center = T)
    # y <- y - ym
  } else {
    ZZ <- t(Z)%*%Z
    SZy <- solve(ZZ,t(Z)%*%y)
    SZX <- solve(ZZ,t(Z)%*%X)
    # X <- X - Z%*%SZX
    # y <- y - Z%*%SZy
  }


  if(n>=p){
    XX <- t(X)%*%X
  } else{
    XX <- X%*%t(X)
  }

  eigenXX <- eigen(XX)
  eVal <- eigenXX$values
  eVec <- eigenXX$vectors

  #initialize
  if(is.null(se2)) {se2 <- drop(var(y))}
  if(is.null(sb2)) {sb2 <- drop(var(y))}
  mu <- matrix(0,p,1)
  beta0 <- SZy - SZX %*% mu
  lb <- rep(0,maxIter)
  lb[1] <- -Inf


  y_bar <- y - Z %*% beta0
  gam <- 1
  for(iter in 2:maxIter){
    #E-step
    # S <- 1/se2 * XX + diag(1/sb2,p)
    D <- eVal/se2 + 1/sb2
    if(n>=p){
      mu <- 1/se2 * eVec %*% (t(eVec)%*%(t(X)%*%y_bar) / D)
    } else {
      mu <- 1/se2 * t(X) %*% (eVec %*% (t(eVec)%*%y_bar / D))
    }

    # mu <- 1/se2 * eVec %*% (t(eVec)%*%Xy / D)
    Xmu <- X %*% mu
    y_Xmu2 <- sum((y_bar-Xmu)^2)

    #if n<d, add zero eigenvalues
    if(n>=p){
      D1 <- D
    } else {
      D1 <- c(D,rep(1/sb2,p-n))
    }

    #evaluate lower bound
    lb[iter] <- get_lb_PXEM(n,sb2,se2,D1,mu,y_Xmu2)
    # loglik <- dmvnorm(t(y_bar),mean=rep(0,n),sigma=X%*%t(X)*sb2+diag(n)*se2,log=TRUE)

    if(verbose){
      cat(iter,"-th iteration, lower bound = ",lb[iter]," ,diff=",lb[iter]-lb[iter-1],",sb2=",sb2,",se2=",se2,"\n")
    }
    if(abs(lb[iter]-lb[iter-1])<tol){
      lb <- lb[1:iter]
      break
    }

    #M-step
    gam <- sum(y_bar*Xmu) / (sum(Xmu^2) + sum(eVal/D))

    beta0 <- SZy - SZX %*% mu * gam
    y_bar <- y - Z %*% beta0

    se2 <- sum((y_bar-Xmu*gam)^2)/n + gam^2 * sum(eVal/D)/n
    sb2 <- sum((mu)^2)/p + sum(1/D1)/p
    # safe guard for small variance components
    sb2 <- ifelse(sb2<1e-6,1e-6,sb2)
    se2 <- ifelse(se2<1e-6,1e-6,se2)

    #Reduciton-step
    sb2 <- gam^2 * sb2
    gam <- 1
  }

  w  <- c(beta0,mu)
  gamma <- p - sum(1/D1)/sb2

  invSigy <- solve(sb2*K+se2*diag(n))
  invSigyK <- invSigy%*%K
  FIM <- matrix(0,2,2)
  FIM[1,1] <- sum(invSigyK^2) / 2
  FIM[2,2] <- sum(invSigy^2) / 2
  FIM[1,2] <- FIM[2,1] <- sum(invSigyK*invSigy) / 2
  covSig <- solve(FIM)  #inverse of FIM

  bayesReg <- list(beta0=beta0,sb2=sb2,se2=se2,mu=mu,gamma=gamma,iter=iter,covSig=covSig,lb=lb)
  attr(bayesReg,"class") <- "VCM_PXEM"
  bayesReg

}



####################################################################################################
get_lb_PXEM <- function(n,sb2,se2,D,mu,y_Xmu2){
  p <- length(mu)
  E <- y_Xmu2/(2*se2) + sum(mu^2)/(2*sb2)
  ret <- - p*log(sb2)/2 - n*log(se2)/2 - E - sum(log(D))/2 - n/2*log(2*pi)
  drop(ret)
}
