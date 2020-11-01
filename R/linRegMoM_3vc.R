linRegMoM_3vc <- function(K1,K2,y,Z=NULL){
  n <- length(y)

  if(is.null(Z)){
    Z <- matrix(1,n,1)
    M <- diag(n) - matrix(1/n,n,n)
    y <- y - mean(y)
    MK1 <- K1
    MK2 <- K2
  } else{
    Z <- cbind(1,Z)
    M <- diag(n) - Z%*%solve(t(Z)%*%Z)%*%t(Z)
    y <- M %*%y
    MK1 <- M%*%K1
    MK2 <- M%*%K2
  }
  q <- ncol(Z)

  trK1 <- sum(diag(MK1))
  trK2 <- sum(diag(MK2))
  trK1K2 <- sum(MK1*MK2)
  trK12 <- sum(MK1^2)
  trK22 <- sum(MK2^2)

  S <- matrix(c(trK12, trK1K2, trK1, trK1K2, trK22, trK2, trK1, trK2,n-q),3,3)
  c <- c(t(y)%*%MK1%*%y, t(y)%*%MK2%*%y, sum(y^2))

  invS <- solve(S)
  sigma <- invS %*% c

  covB <- matrix(0,3,3)
  if(q==1){
    Sigma <- sigma[1]*K1 + sigma[2]*K2 + sigma[3]*M

    K1S <- K1%*%Sigma
    K2S <- K2%*%Sigma

    covB[1,1] <- sum(K1S^2) * 2
    covB[2,2] <- sum(K2S^2) * 2
    covB[3,3] <- sum(Sigma^2) * 2
    covB[1,2] <- covB[2,1] <- sum(K1S*K2S) * 2
    covB[1,3] <- covB[3,1] <- sum(K1S*Sigma) * 2
    covB[2,3] <- covB[3,2] <- sum(K2S*Sigma) * 2
  } else{
    MS <- sigma[1]*MK1%*%M + sigma[2]*MK2%*%M + sigma[3]*M
    MK1MS <- MK1 %*% MS
    MK2MS <- MK2 %*% MS


    covB[1,1] <- sum(MK1MS^2) * 2
    covB[2,2] <- sum(MK2MS^2) * 2
    covB[3,3] <- sum(MS^2) * 2
    covB[1,2] <- covB[2,1] <- sum(MK1MS*MK2MS) * 2
    covB[1,3] <- covB[3,1] <- sum(MK1MS*MS) * 2
    covB[2,3] <- covB[3,2] <- sum(MK2MS*MS) * 2
  }



  # Exact test using Davies method
  covSig <- invS %*% covB %*% invS

  s12 <- sigma[1]
  s22 <- sigma[2]
  se2 <- sigma[3]

  PVE <- matrix(NA,1,3)
  colnames(PVE) <- c("PVE1","PVE2","PVE")

  var_total <- s12*trK1 + s22*trK2 + se2*(n-q)
  PVE[1,] <- c(s12*trK1/var_total, s22*trK2/var_total, (s12*trK1+s22*trK2)/var_total)

  # compute H
  H <- invS[1,1]*MK1%*%M+ invS[1,2]*MK2%*%M + invS[1,3]*M

  ret <- list(s12=s12,s22=s22,se2=se2,K1=K1,K2=K2,PVE=PVE)

}
