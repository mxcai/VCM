# Evaluate the performance of bayes regression using cv
cvperf.linReg <- function(X, y, Z=NULL, nfolds=10, alg="PXEM", verbose=F, ...){
  p <- ncol(X)
  n <- nrow(X)
  # decide the cv assignments
  if(nfolds == n) {
    idx <- sample(1:n)
  } else {
    idx <- ceiling(sample(1:n)/n*nfolds)
  }

  testErr  <- rep(0,nfolds)
  R2       <- rep(0,nfolds)

  # report settings
  message("Info: Algorithm used: ", alg)
  message("Info: Number of variables: ", p)
  message("Info: Sample size: ", n)
  message("Info: Number of cv folds: ", nfolds)

  cat("start cv process......... total",nfolds,"validation sets \n")

  for(i in 1:nfolds) {
    cat(i,"-th validation set... \n")

    X_train <- X[idx!=i,]
    y_train <- y[idx!=i]
    Z_train <- Z[idx!=i,]

    X_test  <- matrix(X[idx==i,],ncol = p)
    y_test  <- y[idx==i]
    Z_test  <- Z[idx==i,]
    if(!is.null(Z)) Z_test <- matrix(Z_test,ncol = ncol(Z))

    if(alg == "MoM"){
      fit <- linRegMoM(X_train,y_train,Z_train,verbose=verbose,...)
    } else if(alg == "EM"){
      fit <- linRegMM(X_train,y_train,Z_train,verbose=verbose,...)
    } else if(alg == "PXEM"){
      fit <- linRegPXEM(X_train,y_train,Z_train,verbose=verbose,...)
    }


    testErr[i]  <- mean((y_test-predict(fit,X_test,Z_test))^2)
    R2[i]       <- 1 - testErr[i] / mean((y_test-mean(y_train))^2)
  }
  cvm  <- mean(testErr)
  cvsd <- sd(testErr)/sqrt(nfolds)
  return(list(testErr=testErr,cvm=cvm,cvsd=cvsd,R2=R2))
}
