library(ggplot2)
library(reshape2)
n <- 1000
d <- 1000
q <- 5

sg2 <- 0.1
se2 <- 1

X <- matrix(rnorm(n*d),n,d)
X <- scale(X)/sqrt(d)
w  <- c(rnorm(d,0,sqrt(sg2)))
y0 <- X%*%w

Z <- matrix(rnorm(n*q),n,q)
omega <- c(1:5)
y1 <- Z%*%omega

y  <- y0 + y1 + sqrt(se2)*rnorm(n)

# fit PXEM
fit_PXEM <- linRegPXEM(X=X,y=y,Z = Z,tol = 1e-6,maxIter =200)
yhat_PX <- predict(fit_PXEM,X,Z)
mean((yhat_PX-y)^2)

# fit MM
fit_MM <- linRegMM(X=X,y=y,Z=Z,tol=1e-6,maxIter = 200)
yhat_MM <- predict(fit_MM,X,Z)
mean((yhat_MM-y)^2)

# fit MoM
fit_MoM <- linReg_MoM(X=X,Z=Z,y=y)

lb <- cbind(PXEM=fit_PXEM$lb[-1],MM=fit_MM$lb[-1])

lb <- melt(lb)
names(lb) <- c('iter','method','likelihood')
P <- ggplot(lb,aes(x=iter,y=likelihood,color=method)) + geom_line()
P

