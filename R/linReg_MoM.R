# VCM solved by Method of Moments
linReg_MoM <- function(X,y,Z=NULL){
  p <- ncol(X)
  n <- length(y)

  ym <- mean(y)
  X <- scale(X)/sqrt(p)
  K <- X %*% t(X)

  if(is.null(Z)){
    Z <- matrix(1,n,1)
    M <- diag(n) - matrix(1/n,n,n)
    y <- y - mean(y)
    MK <- K
  } else{
    Z <- cbind(1,Z)
    M <- diag(n) - Z%*%solve(t(Z)%*%Z)%*%t(Z)
    y <- M %*%y
    MK <- M%*%K
  }
  q <- ncol(Z)

  trK <- sum(diag(MK))
  trK2 <- sum(MK^2)

  S <- matrix(c(trK2, trK, trK, n-q),2,2)
  c <- c(t(y)%*%MK%*%y, sum(y^2))

  invS <- solve(S)
  sigma <- invS %*% c

  covB <- matrix(0,2,2)
  if(q==1){
    Sigma <- sigma[1]*K + sigma[2]*M
    KS <- K%*%Sigma

    covB[1,1] <- sum(KS^2) * 2
    covB[2,2] <- sum(Sigma^2) * 2
    covB[1,2] <- covB[2,1] <- sum(KS*Sigma) * 2
  } else{
    MS <- sigma[1]*MK%*%M + sigma[2]*M
    MKMS <- MK %*% MS

    covB[1,1] <- sum(MKMS^2) * 2
    covB[2,2] <- sum(MS^2) * 2
    covB[1,2] <- covB[2,1] <- sum(MKMS*MS) * 2
  }
  Omega <- sigma[1]*K
  diag(Omega) <- diag(Omega) + sigma[2]

  inv_Omega <- chol2inv(chol(Omega))
  # recover beta0 and mu
  beta0 <- solve(t(Z)%*%inv_Omega%*%Z,t(Z)%*%(inv_Omega%*%y0))

  # Sandwich estimator
  covSig <- invS %*% covB %*% invS
  mu <- sigma[1]*t(X)%*%(inv_Omega%*%(y0-Z%*%beta0))

  mu <- mu/Xsd/sqrt(p)

  sb2 <- sigma[1]
  se2 <- sigma[2]

  var_total <- sb2*trK + se2*(n-q)

  h <- sb2*trK / var_total
  # Delta method
  gh <- c(se2*(n-q)*trK/var_total^2, -sb2*(n-q)*trK/var_total^2)
  se_h <- sqrt(t(gh) %*% covSig %*% gh)

  ret <- list(beta0=beta0, sb2=sb2, se2=se2, mu=mu, K=K, covSig=covSig, h=h, se_h=se_h)

  class(ret) <- c("VCM","MoM")
  ret
}
