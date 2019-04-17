library(ggplot2)
library(reshape2)
n <- 500
d <- 1000

sg2 <- 0.1
se2 <- 1
X <- matrix(rnorm(n*d),n,d)
X <- scale(X)/sqrt(d)

w  <- c(rnorm(d,0,sqrt(sg2)))
y0 <- X%*%w

y  <- y0 + sqrt(se2)*rnorm(n)

fit_PXEM <- linRegPXEM(X=X,y=y,tol = 1e-6,maxIter =2000)
fit_MM <- linRegMM(X=X,y=y,tol=1e-6,maxIter = 2000)
fit_MoM <- linReg_MoM(X=X,y=y)

lb_MM <- fit_MM$lb
lb_PX <- fit_PXEM$lb

lb_MM <- c(lb_MM,rep(NA,ifelse(length(lb_PX)>length(lb_MM),length(lb_PX)-length(lb_MM),0)))
lb_PX <- c(lb_PX,rep(NA,ifelse(length(lb_MM)>length(lb_PX),length(lb_MM)-length(lb_PX),0)))

# maxiter <- min(c(length(fit_MM$lb),length(fit_PXEM$lb)))-1
lb <- cbind(PXEM=lb_PX[2:length(lb_PX)],MM=lb_MM[2:length(lb_MM)])

lb <- melt(lb)
names(lb) <- c('iter','method','likelihood')
P <- ggplot(lb,aes(x=iter,y=likelihood,color=method)) + geom_line()
P
