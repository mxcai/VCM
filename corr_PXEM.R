# SNP-gene-phenotype prediction
corr_PXEM <- function(y,z,X1,X2,U1=NULL,U2=NULL,maxIter=1500,tol=1e-8){
  p <- ncol(X1)
  n1 <- length(y)
  n2 <- length(z)

  if(is.null(U1)){
    U1 <- matrix(1,n1,1)
  } else{
    U1 <- cbind(1,U1)
  }
  if(is.null(U2)){
    U2 <- matrix(1,n2,1)
  } else {
    U2 <- cbind(1,U2)
  }

  #means of X, y and z
  X1m <- colMeans(X1)
  X2m <- colMeans(X2)
  ym <- mean(y)
  zm <- mean(z)


  #calculate sth in advance
  if(ncol(U1)==1){
    SU1y <- ym
    SU1X <- X1m
  } else {
    UU1 <- t(U1)%*%U1
    SU1y <- solve(UU1,t(U1)%*%y)
    SU1X <- solve(UU1,t(U1)%*%X1)
  }
  if(ncol(U2)==1){
    SU2z <- zm
    SU2X <- X2m
  } else {
    UU2 <- t(U2)%*%U2
    SU2z <- solve(UU2,t(U2)%*%z)
    SU2X <- solve(UU2,t(U2)%*%X2)
  }

  XX1 <- t(X1)%*%X1
  XX2 <- t(X2)%*%X2
  # X1y <- t(X1)%*%y
  # X2z <- t(X2)%*%z

  eigenXX1 <- eigen(XX1)
  eVal1 <- eigenXX1$values
  eVec1 <- eigenXX1$vectors

  eigenXX2 <- eigen(XX2)
  eVal2 <- eigenXX2$values
  eVec2 <- eigenXX2$vectors

  #initialize
  sy2 <- drop(var(y))
  s1 <- sy2 / p
  sz2 <- drop(var(z))
  s2 <- sz2 / p

  mu1 <- matrix(0,p,1)
  w1 <- SU1y - SU1X %*% mu1
  mu2 <- matrix(0,p,1)
  w2 <- SU2z - SU2X %*% mu2

  lb0 <- rep(0,maxIter)
  lb0[1] <- -Inf

  #Fit Null model as warm-start: set rho = 0
  #evaluate LSE of alpha0 and sz2
  y_bar <- y - U1 %*% w1
  z_bar <- z - U2 %*% w2

  gamma1 <- 1
  gamma2 <- 1
  for(iter in 2:maxIter){
    #E-step
    # invS <- 1/sy2 * XX1 + diag(1/sb2,p)
    D1 <- eVal1/sy2 + 1/s1
    D2 <- eVal2/sz2 + 1/s2
    mu1 <- 1/sy2 * eVec1 %*% (t(eVec1)%*%(t(X1)%*%y_bar) / D1)
    mu2 <- 1/sz2 * eVec2 %*% (t(eVec2)%*%(t(X2)%*%z_bar) / D2)

    X1mu <- X1%*%mu1
    X2mu <- X2%*%mu2

    y_Xmu2 <- sum((y_bar-X1mu)^2)
    z_Xmu2 <- sum((z_bar-X2mu)^2)

    #evaluate lower bound
    # lb0[iter] <- get_lb1_bv_full(n1,n2,y_bar,z_bar,X1%*%t(X1),X2%*%t(X2),X1%*%t(X2),diag(c(s1,s2)),sy2,sz2)
    lb0[iter] <- get_lb0_bv(n1,n2,s1,s2,sy2,sz2,D1,D2,mu1,mu2,y_Xmu2,z_Xmu2)
    cat("Null model at ",iter,"-th iteration, lower bound = ",lb0[iter]," ,diff=",lb0[iter]-lb0[iter-1],",s1=",s1,",s2=",s2,",sy2=",sy2,",sz2=",sz2,"\n")
    if(abs(lb0[iter]-lb0[iter-1])<1e-8*abs(lb0[iter])){
      lb0 <- lb0[1:iter]
      break
    }

    #M-step
    gamma1 <- sum(y_bar*X1mu) / (sum(X1mu^2) + sum(eVal1/D1))
    gamma2 <- sum(z_bar*X2mu) / (sum(X2mu^2) + sum(eVal2/D2))

    w1 <- SU1y - SU1X %*% mu1 * gamma1
    y_bar <- y - U1 %*% w1

    w2 <- SU2z - SU2X %*% mu2 * gamma2
    z_bar <- z - U2 %*% w2

    sy2 <- sum((y_bar-X1mu*gamma1)^2)/n1 + gamma1^2 * sum(eVal1/D1)/n1
    sz2 <- sum((z_bar-X2mu*gamma2)^2)/n2 + gamma2^2 * sum(eVal2/D2)/n2

    s1 <- sum((mu1)^2)/p + sum(1/D1)/p
    s2 <- sum((mu2)^2)/p + sum(1/D2)/p

    #Reduciton-step
    s1 <- gamma1^2 * s1
    s2 <- gamma2^2 * s2

    gamma1 <- 1
    gamma2 <- 1
  }
  w1 <- SU1y - SU1X %*% mu1
  w2 <- SU2z - SU2X %*% mu2
  # return(list(w1=w1,w2=w2,mu1=mu1,mu2=mu2,sy2=sy2,sz2=sz2,s1=s1,s2=s2,lb=lb0))

  #Fit Full model:
  lb1 <- rep(0,maxIter)
  lb1[1] <- tail(lb0,1)

  y_bar <- y - U1 %*% w1
  z_bar <- z - U2 %*% w2
  Sigma <- diag(c(s1,s2))
  Sj <- matrix(0,2,2)

  XX <- matrix(0,2*p,2*p)

  Gamma <- diag(2)
  for(iter in 2:maxIter){
    #E-step
    XX[1:p,1:p] <- XX1/sy2
    XX[(p+1):(2*p),(p+1):(2*p)] <- XX2/sz2

    cSigma <- chol(Sigma)
    invS <- XX + (chol2inv(cSigma) %x% diag(p))
    cinvS <- chol(invS)
    S <- chol2inv(cinvS)
    Sj[1,1] <- sum(diag(S)[1:p])
    Sj[2,2] <- sum(diag(S)[(p+1):(2*p)])
    Sj[1,2] <- Sj[2,1] <- sum(diag(S[1:p,(p+1):(2*p)]))

    mu <- S %*% c(1/sy2 * t(X1)%*%y_bar, 1/sz2 * t(X2)%*%z_bar)
    muj <- matrix(mu,p,2)

    X1mu1 <- X1%*%muj[,1]
    X1mu2 <- X1%*%muj[,2]
    X2mu1 <- X2%*%muj[,1]
    X2mu2 <- X2%*%muj[,2]

    trXX1S1 <- sum(diag(XX1%*%S[1:p,1:p]))
    trXX1S2 <- sum(diag(XX1%*%S[(p+1):(2*p),(p+1):(2*p)]))
    trXX1S12 <- sum(diag(XX1%*%S[1:p,(p+1):(2*p)]))

    trXX2S1 <- sum(diag(XX2%*%S[1:p,1:p]))
    trXX2S2 <- sum(diag(XX2%*%S[(p+1):(2*p),(p+1):(2*p)]))
    trXX2S12 <- sum(diag(XX2%*%S[1:p,(p+1):(2*p)]))

    y_Xmu2 <- sum((y_bar-X1mu1)^2)
    z_Xmu2 <- sum((z_bar-X2mu2)^2)



    #evaluate lower bound
    # lb1[iter] <- get_lb1_bv_full(n1,n2,y_bar,z_bar,X1%*%t(X1),X2%*%t(X2),X1%*%t(X2),Sigma,sy2,sz2)
    lb1[iter] <- get_lb1_bv(n1,n2,cSigma,sy2,sz2,cinvS,muj,y_Xmu2,z_Xmu2)
    cat("Full model at ",iter,"-th iteration, lower bound = ",lb1[iter]," ,diff=",lb1[iter]-lb1[iter-1],",s1=",Sigma[1,1],",s2=",Sigma[2,2],",rho=",Sigma[1,2]/sqrt(Sigma[1,1])/sqrt(Sigma[2,2]),",sy2=",sy2,",sz2=",sz2,"\n")
    if(iter>2 & abs(lb1[iter]-lb1[iter-1])<tol*abs(lb1[iter])){
      lb1 <- lb1[1:iter]
      break
    }

    #M-step
    H1 <- H2 <- matrix(0,2,2)
    f1 <- f2 <- matrix(0,2,1)

    H1[1,1] <- sum(X1mu1^2) + trXX1S1
    H1[1,2] <- H1[2,1] <- sum(X1mu2*X1mu1) + trXX1S12
    H1[2,2] <- sum(X1mu2^2) + trXX1S2
    f1[1] <- sum(X1mu1*y_bar)
    f1[2] <- sum(X1mu2*y_bar)

    H2[1,1] <- H2[2,2] <- sum(X2mu1*X2mu2) + trXX2S12
    H2[1,2] <- sum(X2mu2^2) + trXX2S2
    H2[2,1] <- sum(X2mu1^2) + trXX2S1
    f2[1] <- sum(X2mu2*z_bar)
    f2[2] <- sum(X2mu1*z_bar)

    Gamma_vec <- c(solve(H1,f1),solve(H2,f2))
    Gamma <- matrix(Gamma_vec,2,2,byrow = T)
    # Gamma[1,1] <- sum(y_bar*X1mu1) / (sum(X1mu1^2) + trXX1S1)
    # Gamma[2,2] <- sum(z_bar*X2mu2) / (sum(X2mu2^2) + trXX2S2)

    w1 <- SU1y - SU1X %*% (mu1*Gamma[1,1]+mu2*Gamma[1,2])
    y_bar <- y - U1 %*% w1

    w2 <- SU2z - SU2X %*% (mu2*Gamma[2,2]+mu1*Gamma[2,1])
    z_bar <- z - U2 %*% w2

    sy2 <- sum((y_bar-X1mu1*Gamma[1,1]-X1mu2*Gamma[1,2])^2)/n1 + Gamma[1,1]^2*trXX1S1/n1 + 2*Gamma[1,1]*Gamma[1,2]*trXX1S12/n1 + Gamma[1,2]^2*trXX1S2/n1
    sz2 <- sum((z_bar-X2mu2*Gamma[2,2]-X2mu1*Gamma[2,1])^2)/n2 + Gamma[2,2]^2*trXX2S2/n2 + 2*Gamma[2,2]*Gamma[2,1]*trXX2S12/n2 + Gamma[2,1]^2*trXX2S1/n2
    Sigma <- (t(muj)%*%muj+Sj)/p

    #Reduction step
    Sigma <- Gamma %*% Sigma %*% t(Gamma)
    Gamma <- diag(2)
  }

  lr <- 2 * (tail(lb1,1) - tail(lb0,1))

  ret <- list()
  ret$sy2 <- sy2
  ret$sz2 <- sz2
  ret$Sigma <- Sigma

  ret$w1 <- w1
  ret$w2 <- w2

  ret$mu <- mu
  ret$S <- S
  ret$lb0 <- lb0
  ret$lb1 <- lb1

  ret$lr <- lr

  ret
}





####################################################################################################
get_lb0_bv <- function(n1,n2,s1,s2,sy2,sz2,D1,D2,mu1,mu2,y_Xmu2,z_Xmu2){
  p <- length(D1)
  E <- y_Xmu2/(2*sy2) + z_Xmu2/(2*sz2) + sum(mu1^2)/(2*s1) + sum(mu2^2)/(2*s2)
  ret <- - p*log(s1)/2 - p*log(s2)/2 - n1*log(sy2)/2 - n2*log(sz2)/2 - E - sum(log(D1))/2 - sum(log(D2))/2 - (n1/2+n2/2)*log(2*pi)
  drop(ret)
}
get_lb1_bv <- function(n1,n2,cSigma,sy2,sz2,cinvS,muj,y_Xmu2,z_Xmu2){
  p <- nrow(cinvS)/2
  E <- y_Xmu2/(2*sy2) + z_Xmu2/(2*sz2) + sum(diag( t(muj) %*% muj %*% chol2inv(cSigma) ))/2
  ret <- - p*sum(log(diag(cSigma))) - n1*log(sy2)/2 - n2*log(sz2)/2 - E - sum(log(diag(cinvS))) - (n1/2+n2/2)*log(2*pi)
  drop(ret)
}
get_lb1_bv_full <- function(n1,n2,y_bar,z_bar,X1X1,X2X2,X1X2,Sigma,sy2,sz2){
  s1 <- Sigma[1,1]
  s2 <- Sigma[2,2]
  s12 <- Sigma[1,2]
  St <- matrix(0,n1+n2,n1+n2)
  St[1:n1,1:n1] <- X1X1*s1 + sy2*diag(n1)
  St[(n1+1):(n1+n2),(n1+1):(n1+n2)] <- X2X2*s2 + sz2*diag(n2)
  St[1:n1,(n1+1):(n1+n2)] <- X1X2*s12
  St[(n1+1):(n1+n2),1:n1] <- t(St[1:n1,(n1+1):(n1+n2)])

  cSt <- chol(St)
  invcSt <- solve(cSt)
  tmp <- c(y_bar,z_bar)%*%invcSt
  ret <- - sum(log(diag(cSt))) - 0.5 * sum(tmp*tmp) - (n1+n2)*0.5*log(2*pi)
  ret
}
