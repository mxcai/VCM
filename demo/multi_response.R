library(ggplot2)
library(reshape2)
n <- 1000
d <- 1000
q <- 5
nY <- 3

sg2 <- 0.1
se2 <- 1

X <- matrix(rnorm(n*d),n,d)
X <- scale(X)/sqrt(d)
w  <- replicate(nY,rnorm(d,0,sqrt(sg2)))
y0 <- X%*%w

Z <- matrix(rnorm(n*q),n,q)
omega <- matrix(c(11:15,51:55,-2:2),q,nY)
y1 <- Z%*%omega

y  <- y0 + y1 + sqrt(se2)*rnorm(n)

fit_mvPXEM <- mvRegPXEM(X=X,Y=y,Z = Z,tol = 1e-6,maxIter =200)
Yhat_PX <- predict(fit_mvPXEM,X,Z)
colMeans((Yhat_PX)^2)  # MSE for mv-PXEM

fit_PXEM <- linRegPXEM(X=X,y=y[,1],Z = Z,tol = 1e-6,maxIter =200)
Yhat_PX3 <- predict(fit_PXEM,X,Z)
mean((Yhat_PX3)^2)  # MSE for PXEM

fit_mvMM <- mvRegMM(X=X,Y=y,Z = Z,tol = 1e-6,maxIter =200)
Yhat_MM <- predict(fit_mvMM,X,Z)
colMeans((Yhat_MM)^2)  # MSE for mv-MM

fit_MM <- linRegMM(X=X,y=y[,3],Z=Z,tol=1e-6,maxIter = 200)
Yhat_MM3 <- predict(fit_MM,X,Z)
mean((Yhat_MM3)^2)  # MSE for MM

