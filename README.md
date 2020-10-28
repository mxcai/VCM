VCM
===========
R package 'VCM' contains three approaches for solving variance components model: PX-EM algorithm, MM algorithm and Method of Moments.

Installation
===========

To install VCM, use the following code in R console:

```{r}
#install.packages("devtools")
library(devtools)
install_github("mxcai/VCM")
```

Usage
===========
[The 'VCM' vignette](https://github.com/mxcai/VCM/blob/master/inst/doc/vignette.pdf) provides details of derivation and a quick start for the usage of the package.

Variance component estimation under model mis-specification
===========
Jiang, J et al. demonstrate that the LMM can produce accurate variance component estimate under mis-specified random effects model. Here we reproduce their results using our VCM package for fitting random effects model and glmnet for fitting lasso model.

We first simulate X with sample size n=500 and p=1000 variables as design matrix. To model the mis-specified setting of the random effects model, we need to generate sparse effects. We consider 3 settings of sparsity level. Specifically, we choose 10, 50, and 100 effects as nonzero by sampling their sizes from normal distributions while fixing the rest effects at zero. For all settings, we set the proportion of variance explained (PVE) at 0.5 (signal-to-noise ratio 1:1). We compare the estimated PVE of noise term.

```{r}
library(glmnet)
library(VCM)
library(ggplot2)
library(RhpcBLASctl)
blas_set_num_threads(8)

n <- 500
p <- 1000

PVEb <- 0.5
PVEe <- 0.5

nonzero <- c(10,50,100)
nrep <- 50

X <- matrix(rnorm(n*p),n,p)
Xs <- scale(X)/sqrt(p)

out <- data.frame()

for(i in 1:length(nonzero)){
  beta <- rep(0,p)

  for (j in 1:nrep){
    beta[1:nonzero[i]] <- rnorm(nonzero[i],0,sqrt(PVEb/nonzero[i]*p))
    y0 <- Xs%*%beta
    y <- y0 + rnorm(n,0,sqrt(PVEe))
    fit_lmm <- linRegMM(X,y)
    fit_lasso <- cv.glmnet(X,y)
    yhat_lasso <- predict(fit_lasso,X,s="lambda.min")
    nz <- sum(coef(fit_lasso,s='lambda.min')!=0)
    PVEe_lasso <- sum((y-yhat_lasso)^2)/(n-nz-1)
    out <- rbind(out,data.frame(nonzero=nonzero[i],method="LMM",PVEe=fit_lmm$se2))
    out <- rbind(out,data.frame(nonzero=nonzero[i],method="LASSO",PVEe=PVEe_lasso))
    cat(i,"-th nonzero, ",j,"-th rep finished.\n")
  }
}
out$nonzero <- as.factor(out$nonzero)

p <- ggplot(out,aes(x=nonzero,y=PVEe,color=method)) + geom_boxplot() + geom_hline(yintercept=PVEe,linetype="dashed")
p
```
![Fig1](https://github.com/mxcai/VCM/blob/master/PVE_LMMvsLASSO.png)

References
==========

Jiang, J., Li, C., Paul, D., Yang, C., & Zhao, H. (2016). On high-dimensional misspecified mixed model analysis in genome-wide association study. The Annals of Statistics, 44(5), 2127-2160.

Liu, C., Rubin, D. B., & Wu, Y. N. (1998). Parameter expansion to accelerate EM: the PX-EM algorithm. Biometrika, 85(4), 755-770.

Foulley, J. L., & Van Dyk, D. A. (2000). The PX-EM algorithm for fast stable fitting of Henderson's mixed model. Genetics Selection Evolution, 32(2), 143.

Zhou, H., Hu, L., Zhou, J., & Lange, K. (2018). MM algorithms for variance components models. Journal of Computational and Graphical Statistics, (just-accepted), 1-30.

Wu, Y., & Sankararaman, S. (2018). A scalable estimator of SNP heritability for biobank-scale data. Bioinformatics, 34(13), i187-i194.


Development
==========

This R package is developed by Mingxuan Cai (mcaiad@ust.hk).


Contact
==========

Please contact Mingxuan Cai (mcaiad@ust.hk) or Prof. Can Yang (macyang@ust.hk) if any enquiry.
