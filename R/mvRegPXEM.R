# ridge regression implemented by PXEM
mvRegPXEM <- function(X,Y,Z=NULL,maxIter=1500,tol=1e-6,se2_init=NULL,sb2_init=NULL,verbose=T){
  Y <- as.matrix(Y)
  nY <- ncol(Y)
  p <- ncol(X)
  n <- nrow(Y)

  Ym <- colMeans(Y)
  Y <- scale(Y,center = T,scale = F)  # centerize y
  Xm <- colMeans(X)
  X <- scale(X,center = T,scale = F)  # centerize X
  Xsd <- sqrt(colMeans(X^2))
  X <- t(t(X)/Xsd)/sqrt(p)  # scale X

  if(is.null(Z)){
    Z <- matrix(1,n,1)
  } else{
    Z <- cbind(1,Z)
  }

  #calculate sth in advance
  if(ncol(Z)==1){
    SZY <- rep(0,nY)
    SZX <- rep(0,p)
    # X <- scale(X,scale = F,center = T)
    # y <- y - ym
  } else {
    ZZ <- t(Z)%*%Z
    SZY <- solve(ZZ,t(Z)%*%Y)
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

  beta0_all <- matrix(0,ncol(Z),nY)
  sb2_all <- rep(0,nY)
  se2_all <- rep(0,nY)
  mu_all <- matrix(0,p,nY)
  gamma_all <- rep(0,nY)
  iter_all  <- rep(0,nY)
  covSig_all <- vector("list",nY)
  lb_all <- vector("list",nY)

  for(k in 1:nY){
    y <- Y[,k]
    ym <- Ym[k]
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

      #Reduciton-step
      sb2 <- gam^2 * sb2
      gam <- 1
      # safe guard for small variance components
      sb2 <- ifelse(sb2<1e-6,1e-6,sb2)
      se2 <- ifelse(se2<1e-6,1e-6,se2)
    }

    gamma <- p - sum(1/D1)/sb2

    # recover beta0 and mu
    mu <- mu/Xsd/sqrt(p)
    if(ncol(Z)==1){
      beta0 <- beta0 + ym - colSums(mu*Xm)
    } else{
      beta0 <- beta0 + solve(ZZ,colSums(Z)*(ym-sum(mu*Xm)))  # = beta0 - (Z^TZ)^-1Z^T[(y^bar-X^bar%*%mu),...,(y^bar-X^bar%*%mu)]
    }

    FIM <- matrix(0,2,2)

    # When n>=p, we use the dual form: X^TOmega^-1 = Lam^-1X^T, where Omega^-1=(sb2*XX^T + se2*I_n)^-1, Lam^-1=(sb2*X^TX + se2*I_p)^-1
    if(n>=p){                                        # work on dimension p if n>=p
      invSigy_p <- eVec%*%(1/(eVal*sb2+se2)*t(eVec))        # p by p dual form: Lam^-1=(se2*I_p + sb2*X^TX)^-1
      invSigyK_p <- invSigy_p%*%XX                          # Lam^-1(X^TX)
      FIM[1,1] <- sum(invSigyK_p^2) / 2                     # tr(Omega^-1XX^TOmega^-1XX^T)=tr(Lam^-1X^TXLam^-1X^TX)
      FIM[2,2] <- n/se2^2 -                                 # dual form of tr(Omega^-2)
        2*sum(diag(invSigyK_p))*sb2/se2^2 +
        sum(invSigyK_p^2)*sb2^2/se2^2
      FIM[1,2] <- FIM[2,1] <- sum(invSigyK_p*invSigy_p) / 2 # tr(Omega^-1XX^TOmega^-1) = tr(Lam^-2X^TX)
    } else {                                         # work on dimension n if n<p
      invSigy <- eVec%*%(1/(eVal*sb2+se2)*t(eVec))          #  n by n original form Omega^-1=(se2*I_n + sb2*XX^T)^-1 <=> solve(sb2*K+se2*diag(n))
      invSigyK <- invSigy%*%XX                              # Omega^-1(XX^T)
      FIM[1,1] <- sum(invSigyK^2) / 2                       # tr(Omega^-1XX^TOmega^-1XX^T)
      FIM[2,2] <- sum(invSigy^2) / 2                        # tr(Omega^-2)
      FIM[1,2] <- FIM[2,1] <- sum(invSigyK*invSigy) / 2     # tr(Omega^-1XX^TOmega^-1)
    }
    # invSigy <- solve(sb2*XX+se2*diag(n))
    # invSigyK <- invSigy%*%XX
    # FIM[1,1] <- sum(invSigyK^2) / 2
    # FIM[2,2] <- sum(invSigy^2) / 2
    # FIM[1,2] <- FIM[2,1] <- sum(invSigyK*invSigy) / 2
    covSig <- solve(FIM)  #inverse of FIM

    beta0_all[,k] <- beta0
    sb2_all[k] <- sb2
    se2_all[k] <- se2
    mu_all[,k] <- mu
    gamma_all[k] <- gamma
    iter_all[k]  <- iter
    covSig_all[[k]] <- covSig
    lb_all[[k]] <- lb
  }


  bayesReg <- list(beta0=beta0_all,sb2=sb2_all,se2=se2_all,mu=mu_all,gamma=gamma_all,iter=iter_all,covSig=covSig_all,lb=lb_all)
  class(bayesReg) <- c("mvVCM","PXEM")
  bayesReg

}



####################################################################################################
get_lb_PXEM <- function(n,sb2,se2,D,mu,y_Xmu2){
  p <- length(mu)
  E <- y_Xmu2/(2*se2) + sum(mu^2)/(2*sb2)
  ret <- - p*log(sb2)/2 - n*log(se2)/2 - E - sum(log(D))/2 - n/2*log(2*pi)
  drop(ret)
}
