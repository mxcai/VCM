library(ggplot2)
library(reshape2)
n <- 1000
d <- 1000

sg2 <- 0.1
se2 <- 1
X <- matrix(rnorm(n*d),n,d)
X <- scale(X)/sqrt(d)

w  <- c(rnorm(d,0,sqrt(sg2)))
y0 <- X%*%w

y  <- y0 + sqrt(se2)*rnorm(n)

fit_PXEM <- linRegPXEM(X=X,y=y,tol = 1e-6,maxIter =200)
fit_MM <- linRegMM(X=X,y=y,tol=1e-6,maxIter = 200)
fit_MoM <- linReg_MoM(X=X,y=y)

lb <- cbind(PXEM=fit_PXEM$lb[-1],MM=fit_MM$lb[-1])

lb <- melt(lb)
names(lb) <- c('iter','method','likelihood')
P <- ggplot(lb,aes(x=iter,y=likelihood,color=method)) + geom_line()
P

