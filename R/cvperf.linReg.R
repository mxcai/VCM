# Evaluate the performance of bayes regression using cv
cvperf.linReg <- function(X=NULL, y, Z=NULL, nfolds=10, alg="PXEM", verbose=F, ...){
  if(is.null(Z)&is.null(X)){
    stop("Either X or Z should be provided.")
  }

  p <- ifelse(is.null(X),0,ncol(X))
  n <- length(y)
  q <- ifelse(is.null(Z),1,ncol(Z)+1)
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
  message("Info: Number of variables with random effects: ", p)
  message("Info: Number of variables with fixed effects (include intercept): ", q)
  message("Info: Sample size: ", n)
  message("Info: Number of cv folds: ", nfolds)

  cat("start cv process......... total",nfolds,"validation sets \n")

  for(i in 1:nfolds) {
    cat(i,"-th validation set... \n")

    y_train <- y[idx!=i]
    y_test  <- y[idx==i]

    Z_train <- Z[idx!=i,]
    Z_test  <- Z[idx==i,]
    if(!is.null(Z)) {
      Z_test <- matrix(Z_test,ncol = ncol(Z))
      colnames(Z_test) <- colnames(Z_train)
    }

    if(is.null(X)){
      # fixed effects model
      fit <- lm(y~.,data = data.frame(y=y_train,Z=Z_train))

      testErr[i]  <- mean((y_test-predict(fit,data.frame(Z=Z_test)))^2)
    } else{
      # mixed effects model
      X_train <- X[idx!=i,]

      X_test  <- matrix(X[idx==i,],ncol = p)

      if(alg == "MoM"){
        fit <- linRegMoM(X_train,y_train,Z_train,verbose=verbose,...)
      } else if(alg == "EM"){
        fit <- linRegMM(X_train,y_train,Z_train,verbose=verbose,...)
      } else if(alg == "PXEM"){
        fit <- linRegPXEM(X_train,y_train,Z_train,verbose=verbose,...)
      }
      testErr[i]  <- mean((y_test-predict(fit,X_test,Z_test))^2)
    }
    R2[i]       <- 1 - testErr[i] / mean((y_test-mean(y_train))^2)
  }
  cvm  <- mean(testErr)
  cvsd <- sd(testErr)/sqrt(nfolds)
  return(list(testErr=testErr,cvm=cvm,cvsd=cvsd,R2=R2))
}
