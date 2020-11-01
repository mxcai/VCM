linRegMM_3vc <- function(K1,K2,y,Z=NULL,maxIter=1500,tol=1e-8,verbose=T,s12=NULL,s22=NULL,se2=NULL){
  n <- length(y)

  if(is.null(Z)){
    Z <- matrix(1,n,1)
  } else{
    Z <- cbind(1,Z)
  }

  #means of X, y and z
  ym <- mean(y)


  #calculate sth in advance
  if(ncol(Z)==1){
    SZy <- ym
  } else {
    ZZ <- t(Z)%*%Z
    SZy <- solve(ZZ,t(Z)%*%y)
  }

  # eigenK <- eigen(K)
  # eVal <- eigenK$values
  # eVec <- eigenK$vectors
  #
  # U <- eVec
  #
  # yt <- t(U)%*%y
  # Zt <- t(U)%*%Z

  #initialize
  if(is.null(se2)) se2 <- var(drop(y))
  if(is.null(s12)) s12 <- var(drop(y)) #/ 90
  if(is.null(s22)) s22 <- var(drop(y)) #/ 90
  # mu <- matrix(0,p,1)
  # beta0 <- SZy - SZX %*% mu
  lb0 <- rep(0,maxIter)
  lb0[1] <- -Inf

  for(iter in 2:maxIter){

    # d <- 1/(sb2*eVal + sy2)
    # beta0 <- solve(t(Zt)%*%diag(d)%*%Zt) %*% t(Zt) %*% diag(d) %*% yt
    # res <- yt - Zt %*% beta0
    #
    # sb2 <- sb2 * sqrt(sum(res^2 * d^2 * eVal) / sum(d * eVal))
    # sy2 <- sy2 * sqrt(sum(res^2 * d^2) / sum(d))

    Omega <- s12*K1 + s22*K2 + se2*diag(n)
    chol_Omega <- chol(Omega)
    inv_Omega <- chol2inv(chol_Omega)

    ZtO <- t(Z)%*%inv_Omega
    beta0 <- solve(ZtO%*%Z) %*% (ZtO%*%y)

    res <- y - Z%*%beta0
    temp <- t(res) %*% inv_Omega

    s12 <- s12 * drop(sqrt(temp%*%K1%*%t(temp) / sum(inv_Omega*K1)))
    s22 <- s22 * drop(sqrt(temp%*%K2%*%t(temp) / sum(inv_Omega*K2)))
    se2 <- se2 * drop(sqrt(sum(temp^2) / sum(diag(inv_Omega))))

    #evaluate lower bound
    lb0[iter] <- get_lb_MM3(chol_Omega,res,temp)

    if(verbose){
      cat(iter,"-th iteration, lower bound = ",lb0[iter]," ,diff=",lb0[iter]-lb0[iter-1],",s12=",s12,",s22=",s22,",se2=",se2,"\n")
    }
    if(abs(lb0[iter]-lb0[iter-1])<(tol*abs(lb0[iter]))){
      lb0 <- lb0[1:iter]
      break
    }
  }

  # calculate covariance of parameter estimates by inverse FIM
  invOmegaK1 <- inv_Omega %*% K1
  invOmegaK2 <- inv_Omega %*% K2
  FIM <- matrix(0,3,3)
  FIM[1,1] <- sum(invOmegaK1^2) / 2
  FIM[2,2] <- sum(invOmegaK2^2) / 2
  FIM[3,3] <- sum(inv_Omega^2) / 2
  FIM[1,2] <- FIM[2,1] <- sum(invOmegaK1*invOmegaK2) / 2
  FIM[1,3] <- FIM[3,1] <- sum(invOmegaK1*inv_Omega) / 2
  FIM[2,3] <- FIM[3,2] <- sum(invOmegaK2*inv_Omega) / 2
  covSig <- solve(FIM)  #inverse of FIM

  # calculate standard error of `heritabilities'
  trK1 <- sum(diag(K1))
  trK2 <- sum(diag(K2))
  var_total <- s12*trK1 + s22*trK2 + se2*n

  H <- matrix(0,2,4)
  colnames(H) <- c("PVEg","PVEa","h2","prop")
  rownames(H) <- c("heritability","se")
  H[1,] <- c(s12*trK1/var_total, s22*trK2/var_total, (s12*trK1+s22*trK2)/var_total, s12*trK1/(s12*trK1 + s22*trK2))

  # calculate h_R = (s12*trK1) / (var_total)
  gh <- c((s22*trK1*trK2+se2*n*trK1)/var_total^2, -s12*trK1*trK2/var_total^2, -s12*n*trK1/var_total^2)
  H[2,1] <- sqrt(t(gh) %*% covSig %*% gh)

  # calculate h_D = (s22*trK2) / (var_total)
  gh <- c(-s22*trK1*trK2/var_total^2, (s12*trK1*trK2+se2*n*trK2)/var_total^2, -s22*trK2*n/var_total^2)
  H[2,2] <- sqrt(t(gh) %*% covSig %*% gh)

  # calculate h_all = (s12*trK1 + s22*trK2) / (var_total)
  gh <- c(se2*n*trK1/var_total^2, se2*n*trK2/var_total^2, -(s12*trK1*n+s22*trK2*n)/var_total^2)
  H[2,3] <- sqrt(t(gh) %*% covSig %*% gh)

  # calculate h_med = s12*trK1 / (s12*trK1 + s22*trK2)
  gh <- c(s22*trK1*trK2/(s12*trK1 + s22*trK2)^2, -s12*trK1*trK2/(s12*trK1 + s22*trK2)^2)
  H[2,4] <- sqrt(t(gh) %*% covSig[1:2,1:2] %*% gh)



  bayesRegMM <- list(beta0=beta0,s12=s12,s22=s22,se2=se2,K1=K1,K2=K2,iter=iter,covSig=covSig,PVE=H)

  attr(bayesRegMM,"class") <- "bayesRegMM"
  bayesRegMM
}


####################################################################################################
get_lb_MM3 <- function(chol_Omega,res,temp){
  llh <- -sum(log(diag(chol_Omega))) - temp%*%res/2
}
