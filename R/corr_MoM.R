corr_MoM <- function(y,z,X1,X2){
  n1 <- length(y)
  n2 <- length(z)
  p <- ncol(X1)

  # Here K are defined as p by p since p is smaller
  K1 <- t(X1)%*%X1    # X1'*X1
  K2 <- t(X2)%*%X2    # X2'*X2
  # K12 <- X1%*%t(X2)

  trK1 <- sum(diag(K1))
  trK2 <- sum(diag(K2))

  trK1sq <- sum(K1^2)    # tr(X1*X1') = tr(X1'*X1)
  trK2sq <- sum(K2^2)    # tr(X*X2') = tr(X2'*X2)
  trK12sq <- sum(K1*K2)    # tr(X1*X2'*X2*X1') = tr(X1'*X1*X2'*X2)

  Xy <- t(X1)%*%y
  yKy <- sum((Xy)^2)
  Xz <- t(X2)%*%z
  zKz <- sum((Xz)^2)

  yKz <- sum(Xy*Xz)

  # solve for expression level vc model
  Sy <- matrix(c(trK1sq, trK1, trK1, n1),2,2)
  cy <- c(yKy, sum(y^2))

  invSy <- solve(Sy)
  sigmay <- invSy %*% cy

  # solve for phenotype level vc model
  Sz <- matrix(c(trK2sq, trK2, trK2, n2),2,2)
  cz <- c(zKz, sum(z^2))

  invSz <- solve(Sz)
  sigmaz <- invSz %*% cz

  # estimate covariance
  cov <- yKz / trK12sq

  # estimate correlation coefficient
  rho <- cov/sqrt(sigmay[1]*sigmaz[1])

  # Sandwich estimator for variance of sig_y and sig_1
  Sigmay <- sigmay[1]*K1 + sigmay[2]*diag(p)
  K1S <- K1%*%Sigmay

  covBy <- matrix(0,2,2)
  covBy[1,1] <- sum(K1S^2) * 2
  covBy[2,2] <- sum(Sigmay^2) * 2
  covBy[1,2] <- covBy[2,1] <- sum(K1S*Sigmay) * 2
  covSigy <- invSy %*% covBy %*% invSy

  # Sandwich estimator for variance of sig_z and sig_2
  Sigmaz <- sigmaz[1]*K2 + sigmaz[2]*diag(p)
  K2S <- K2%*%Sigmaz

  covBz <- matrix(0,2,2)
  covBz[1,1] <- sum(K2S^2) * 2
  covBz[2,2] <- sum(Sigmaz^2) * 2
  covBz[1,2] <- covBz[2,1] <- sum(K2S*Sigmaz) * 2
  covSigz <- invSz %*% covBz %*% invSz

  # calculate variance of rho using Delta method
  g <- -0.5*cov*c(1/sqrt(sigmay[1]^3*sigmaz[1]),1/sqrt(sigmay[1]*sigmaz[1]^3))
  var_rho <- sum(g^2*c(covSigy[1,1],covSigz[1,1]))

  Sigma <- matrix(c(sigmay,sigmaz,sqrt(diag(covSigy)),sqrt(diag(covSigz))),2,4,byrow = T)

  ret <- list(Sigma=Sigma, rho=c(rho,sqrt(var_rho)), K1=K1, K2=K2, covSigy=covSigy, covSigz=covSigz)
}
