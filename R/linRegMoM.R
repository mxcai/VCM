# VCM solved by Method of Moments
linRegMoM <- function(X,y,Z=NULL,approx_se=F,rv_approx=F,n_rv=20,seed=100){
  set.seed(seed)
  p <- ncol(X)
  n <- length(y)

  ym <- mean(y)
  Xm <- colMeans(X)
  X <- scale(X,center = T,scale = F)  # centerize X
  Xsd <- sqrt(colMeans(X^2))
  X <- t(t(X)/Xsd)/sqrt(p)  # scale X
  K <- X %*% t(X)
  y0 <- y

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
  trK2 <- ifelse(rv_approx,
                 sum((MK%*%matrix(rnorm(n_rv*n),n,n_rv))^2)/n_rv,
                 sum(MK^2))

  S <- matrix(c(trK2, trK, trK, n-q),2,2)
  c <- c(t(y)%*%MK%*%y, sum(y^2))

  invS <- solve(S)
  sigma <- invS %*% c

  var_rv <- ifelse(rv_approx,
                   sigma[1]^2*trK2/n_rv,
                   0)

  covB <- matrix(0,2,2)
  if(approx_se){
    if(q==1){
      Ky <- K%*%y
      Sy <- sigma[1]*Ky + sigma[2]*y
      SKy <- sigma[1]*K%*%Ky + sigma[2]*Ky

      covB[1,1] <- t(Ky)%*%(SKy) * 2 + var_rv
      covB[2,2] <- t(y)%*%Sy * 2
      covB[1,2] <- covB[2,1] <- t(SKy)%*%y * 2
    } else{
      Ky <- MK%*%y
      Sy <- sigma[1]*Ky + sigma[2]*y
      SKy <- sigma[1]*MK%*%Ky + sigma[2]*Ky

      covB[1,1] <- t(Ky)%*%(SKy) * 2 + var_rv
      covB[2,2] <- t(y)%*%Sy * 2
      covB[1,2] <- covB[2,1] <- t(SKy)%*%y * 2
    }
  } else {
    if(q==1){
      Sigma <- sigma[1]*K + sigma[2]*M
      KS <- K%*%Sigma

      covB[1,1] <- sum(KS^2) * 2 + var_rv
      covB[2,2] <- sum(Sigma^2) * 2
      covB[1,2] <- covB[2,1] <- sum(KS*Sigma) * 2
    } else{
      MS <- sigma[1]*MK%*%M + sigma[2]*M
      MKMS <- MK %*% MS

      covB[1,1] <- sum(MKMS^2) * 2 + var_rv
      covB[2,2] <- sum(MS^2) * 2
      covB[1,2] <- covB[2,1] <- sum(MKMS*MS) * 2
    }
  }

  Omega <- sigma[1]*K
  diag(Omega) <- diag(Omega) + sigma[2]

  inv_Omega <- chol2inv(chol(Omega))
  # recover beta0 and mu
  beta0 <- solve(t(Z)%*%inv_Omega%*%Z,t(Z)%*%(inv_Omega%*%y0))
  mu <- sigma[1]*t(X)%*%(inv_Omega%*%(y0-Z%*%beta0))

  # re-scale beta0 and mu
  mu <- mu/Xsd/sqrt(p)
  beta0[1] <- beta0[1] - colSums(mu*Xm)

  # Sandwich estimator
  covSig <- invS %*% covB %*% invS

  sb2 <- sigma[1]
  se2 <- sigma[2]

  var_total <- sb2*trK + se2*(n-q)

  h <- sb2*trK / var_total
  # Delta method
  gh <- c(se2*(n-q)*trK/var_total^2, -sb2*(n-q)*trK/var_total^2)
  se_h <- sqrt(t(gh) %*% covSig %*% gh)

  ret <- list(beta0=beta0, sb2=sb2, se2=se2, mu=mu, K=K, covSig=covSig, h=h, se_h=se_h)

  rm(.Random.seed, envir=.GlobalEnv)
  class(ret) <- c("VCM","MoM")
  ret
}
