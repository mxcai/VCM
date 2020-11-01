# ridge regression implemented by PXEM
linRegPXEM_3vc <- function(X1,X2,y,Z=NULL,maxIter=1500,tol=1e-6,se2=NULL,s1=NULL,s2=NULL,verbose=T,PX=T){
  p1 <- ncol(X1)
  p2 <- ncol(X2)
  p <- p1+p2
  n <- length(y)

  ym <- mean(y)
  y  <- y-ym  # centerize y

  X1m <- colMeans(X1)
  X1 <- scale(X1,center = T,scale = F)  # centerize X1
  X1sd <- sqrt(colMeans(X1^2))
  X1 <- t(t(X1)/X1sd)/sqrt(p1)  # scale X1

  X2m <- colMeans(X2)
  X2 <- scale(X2,center = T,scale = F)  # centerize X2
  X2sd <- sqrt(colMeans(X2^2))
  X2 <- t(t(X2)/X2sd)/sqrt(p2)  # scale X2

  Xsd <- c(X1sd,X2sd)
  Xm <- c(X1m,X2m)

  X <- cbind(X1,X2)

  # K <- X %*% t(X)

  if(is.null(Z)){
    Z <- matrix(1,n,1)
  } else{
    Z <- cbind(1,Z)
  }

  #calculate sth in advance
  if(ncol(Z)==1){
    SZy <- 0
    SZX <- rep(0,p)
    # X <- scale(X,scale = F,center = T)
    # y <- y - ym
  } else {
    ZZ <- t(Z)%*%Z
    SZy <- solve(ZZ,t(Z)%*%y)
    SZX <- solve(ZZ,t(Z)%*%X)
    # X <- X - Z%*%SZX
    # y <- y - Z%*%SZy
  }


  XX <- t(X)%*%X


  #initialize
  if(is.null(se2)) {se2 <- drop(var(y))/3}
  if(is.null(s1)) {s1 <- drop(var(y))/3}
  if(is.null(s2)) {s2 <- drop(var(y))/3}
  mu <- matrix(0,p,1)
  beta0 <- SZy - SZX %*% mu
  lb <- rep(0,maxIter)
  lb[1] <- -Inf


  y_bar <- y - Z %*% beta0
  gam <- 1
  for(iter in 2:maxIter){
    #E-step

    invSig <- XX/se2
    diag(invSig) <- diag(invSig) + 1/c(rep(s1,p1),rep(s2,p2))
    cinvS <- chol(invSig)
    Sig <- chol2inv(cinvS)

    mu <- Sig %*% (1/se2 * t(X)%*%y_bar)
    Xmu <- X %*% mu
    y_Xmu2 <- sum((y_bar-Xmu)^2)

    #evaluate lower bound
    lb[iter] <- get_lb_PXEM3(n,p1,p2,s1,s2,se2,cinvS,mu,y_Xmu2)
    # lb[iter] <- get_lb_PXEM3(y_bar,X1,X2,s1,s2,se2,n)

    if(verbose){
      cat(iter,"-th iteration, lower bound = ",lb[iter]," ,diff=",lb[iter]-lb[iter-1],",s1=",s1,",s2=",s2,",se2=",se2,"\n")
    }
    if(abs(lb[iter]-lb[iter-1])<tol){
      lb <- lb[1:iter]
      break
    }

    #M-step
    TrXXSig <- sum(XX*Sig)

    if(PX){
      gam <- sum(y_bar*Xmu) / (sum(Xmu^2) + TrXXSig)
    } else{
      gam <- 1
    }


    beta0 <- SZy - SZX %*% mu * gam
    y_bar <- y - Z %*% beta0

    se2 <- sum((y_bar-Xmu*gam)^2)/n + gam^2 * TrXXSig/n
    s1 <- sum((mu[1:p1])^2)/p1 + sum(diag(Sig)[1:p1])/p1
    s2 <- sum((mu[(p1+1):p])^2)/p2 + sum(diag(Sig)[(p1+1):p])/p2

    #Reduciton-step
    s1 <- gam^2 * s1
    s2 <- gam^2 * s2
    gam <- 1
    # safe guard for small variance components
    s1 <- ifelse(s1<1e-6,1e-6,s1)
    s2 <- ifelse(s2<1e-6,1e-6,s2)
    se2 <- ifelse(se2<1e-6,1e-6,se2)
  }

  # recover beta0 and mu
  mu <- mu/Xsd/c(rep(sqrt(p1),p1),rep(sqrt(p2),p2))
  if(ncol(Z)==1){
    beta0 <- beta0 + ym - colSums(mu*Xm)
  } else{
    beta0 <- beta0 + solve(ZZ,colSums(Z)*(ym-sum(mu*Xm)))  # = beta0 - (Z^TZ)^-1Z^T[(y^bar-X^bar%*%mu),...,(y^bar-X^bar%*%mu)]
  }

  bayesReg <- list(beta0=beta0,s1=s1,s2=s2,se2=se2,mu=mu,iter=iter,lb=lb)
  class(bayesReg) <- c("VCM","PXEM")
  bayesReg

}



####################################################################################################
get_lb_PXEM3 <- function(n,p1,p2,s1,s2,se2,cinvS,mu,y_Xmu2){
  p <- length(mu)
  E <- y_Xmu2/(2*se2) + sum(mu^2/c(rep(s1,p1),rep(s2,p2)))/2
  ret <- - p1*log(s1)/2 - p2*log(s2)/2 - n*log(se2)/2 - E - sum(log(diag(cinvS))) - n/2*log(2*pi)
# get_lb_PXEM3 <- function(y_bar,X1,X2,s1,s2,se2,n){
#   ret <- dmvnorm(c(y_bar),mean = rep(0,n),sigma = X1%*%t(X1)*s1+X2%*%t(X2)*s2+diag(se2,n),log = T)
  drop(ret)
}
