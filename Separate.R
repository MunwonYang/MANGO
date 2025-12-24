library(Tlasso)
library(tensr)
library(expm)
library(rTensor)
library(glasso)

source("~/Basic Setting.R")
source("~/helper.R")

Separate<- function(x,val = NULL,nlambda = 20,est.mode = NULL,
                    rate = 1e-3,normalize = T,lambda.list = NULL,
                    thr = 1e-3,maxit = 1e3){
  dimen = dim(x) # dimension of data set
  K = length(dimen) - 1 # order of tensor
  n = dimen[K + 1] # sample size of training set
  n.val <- dim(val)[K+1]
  nvars = prod(dimen) # number of variables
  m.vec = dimen[1:K] # dimension of each observation
  
  if (is.null(est.mode) == TRUE) {
    est.mode = c(1:K)
  }
  
  fit1 <- array(list(),K)
  for(mode in 1:K){
    if(n * prod(m.vec[-mode])<=((m.vec[mode]-1)/2*m.vec[mode])){
      Omega_tilde <- diag(m.vec[mode])
    }else{
      S.array <- array(0,c(m.vec[mode],m.vec[mode],n))
      for(i in 1:n){
        e <- 0
        eval(parse(text = paste("e=x[", paste(rep(",", K), collapse = ""), "i]"))) # assign the ith observation to e
        A <- k_unfold(as.tensor(e),mode)@data
        S.array[,,i] <- A %*% t(A)
      }
      S.mat <- apply(S.array,c(1,2),mean) / prod(m.vec[-mode]);Omega_tilde <- solve(S.mat)
      if (normalize) {
        Omega_tilde= Omega_tilde / norm(Omega_tilde, type = "F")}
    }
    fit1[[mode]] <- Omega_tilde
  }
  
  loglik <- array(list(),K)
  if(is.null(lambda.list)){
    lambda.list <- array(list(),K)
  }
  
  if(!is.null(val)){
    Omega.list = list() # list of \tilde\Omega
    Omega.list.sqrt = list() # list of square root of \tilde\Omega
    for (k in 1:K) {
      Omega.list[[k]] = fit1[[k]]
      Omega.list.sqrt[[k]] = sqrtm(Omega.list[[k]])
    }
    Omega.sqrt.copy = Omega.list.sqrt
    
    for (mode_index in 1:K) {
      k = est.mode[mode_index]
      
      # Calculate \tilde S_k using the training set
      S.array = array(0, c(m.vec[k], m.vec[k], n))
      Omega.list.sqrt[[k]] = diag(m.vec[k]) # set \tilde\Omega_k to identity matrix
      for (i in 1:n) {
        d = 0
        eval(parse(text = paste("d=x[", paste(rep(",", K), collapse = ""), "i]"))) # assign the ith observation to d
        Vi = k_unfold(as.tensor(ttl(as.tensor(d), Omega.list.sqrt,
                                    ms = 1:K
        )@data), m = k)@data
        S.array[, , i] = Vi %*% t(Vi)
      }
      S.mat = apply(S.array, c(1, 2), mean) * m.vec[k] / prod(m.vec) # \tilde S_k
      S.mat1 <- S.mat - diag(diag(S.mat))
      if(is.null(lambda.list[[mode_index]])){
      lambda.list[[mode_index]] <- exp(seq(log(rate * norm(S.mat1,"M")), log(norm(S.mat1,"M")),length.out = nlambda))
      }
      
      # Calculate \tilde S_k using the validation set
      testS.array = array(0, c(m.vec[k], m.vec[k], n.val))
      for (i in 1:n.val) {
        d = 0
        eval(parse(text = paste("d=val[", paste(rep(",", K), collapse = ""), "i]")))
        Vi = k_unfold(as.tensor(ttl(as.tensor(d), Omega.list.sqrt,
                                    ms = 1:K
        )@data), m = k)@data
        testS.array[, , i] = Vi %*% t(Vi)
      }
      testS.mat = apply(testS.array, c(1, 2), mean) * m.vec[k] / prod(m.vec) # \tilde S_k
      Omega.list.sqrt[[k]] = Omega.sqrt.copy[[k]]
      
      # fit model with a sequence of lambdas
      lamk = lambda.list[[mode_index]] # a sequence of candidates for lambda_k
      lam.length = length(lamk)
      loglik2 = rep(0, lam.length)
      for (i in 1:lam.length) {
        Out1 = glasso(S.mat, rho = lamk[i],
                      penalize.diagonal = FALSE, maxit = maxit, thr = thr)
        hat_Omega = Out1$wi;calc <- log(eigen(hat_Omega)$values)
        if(!is.complex(eigen(hat_Omega)$values)){
          loglik2[i] = -tr(testS.mat %*% hat_Omega) + sum(calc)
        }else{
          loglik2[i] <- NA
        }
      }
      loglik[[mode_index]] = loglik2
    }
  }
  lambda.best <- array(NA,K)
  for(mode in 1:K){
    lambda.best[mode] <- lambda.list[[mode]][which.max(loglik[[mode]])]
  }
  
  fit.result <- array(list(),K)
  for(mode_ind in 1:K){
    k = est.mode[mode_ind]
    Omega.list.sqrt = list()
    for (i in 1:K) {
      Omega.list.sqrt[[i]] = sqrtm(fit1[[i]])
    }
    # Calculate \tilde S_k
    S.array = array(0, c(m.vec[k], m.vec[k], n))
    Omega.list.sqrt[[k]] = diag(m.vec[k])
    for (i in 1:n) {
      d = 0
      eval(parse(text = paste("d=x[", paste(rep(",", K), collapse = ""), "i]")))
      Vi = k_unfold(as.tensor(ttl(as.tensor(d), Omega.list.sqrt,
                                  ms = 1:K
      )@data), m = k)@data
      S.array[, , i] = Vi %*% t(Vi)
    }
    S.mat = apply(S.array, c(1, 2), mean) * m.vec[k] / prod(m.vec) # \tilde S_k
    # fit model
    Out1 = glasso(S.mat, rho = lambda.best[mode_ind],
                  penalize.diagonal = FALSE, maxit = maxit, thr = thr)
    hat_Omega = as.matrix(Out1$wi)
    # normalization
    if (normalize) {
      hat_Omega = hat_Omega / norm(hat_Omega, type = "F")
    }
    fit.result[[mode_ind]] = hat_Omega
  }
  
  result <- list(Omegahat = fit.result,lambda.best = lambda.best,
                 loglik = loglik, 
                 lambda.list = lambda.list)
  return(result)
}

cv.Separate <- function(x,nfold = 5,est.mode = NULL,
                        seed = 14,nlambda = c(20,10),rate = 1e-3,
                        normalize = TRUE,thr = 1e-3,maxit = 1e3,
                        control.rate = 1e-2){
  set.seed(seed)
  dimen = dim(x) # dimension of data set
  K = length(dimen) - 1 # order of tensor
  n = dimen[K + 1] # sample size of training set
  nvars = prod(dimen) # number of variables
  m.vec = dimen[1:K] # dimension of each observation
  
  if (is.null(est.mode) == TRUE) {
    est.mode = c(1:K)
  }
  
  fit1 <- array(list(),K)
  for(mode in 1:K){
    if(n * prod(m.vec[-mode])<=((m.vec[mode]-1)/2*m.vec[mode])){
      Omega_tilde <- diag(m.vec[mode])
    }else{
      S.array <- array(0,c(m.vec[mode],m.vec[mode],n))
      for(i in 1:n){
        e <- 0
        eval(parse(text = paste("e=x[", paste(rep(",", K), collapse = ""), "i]"))) # assign the ith observation to e
        A <- k_unfold(as.tensor(e),mode)@data
        S.array[,,i] <- A %*% t(A)
      }
      S.mat <- apply(S.array,c(1,2),mean) / prod(m.vec[-mode]);Omega_tilde <- solve(S.mat)
      if (normalize) {
        Omega_tilde= Omega_tilde / norm(Omega_tilde, type = "F")}
    }
    fit1[[mode]] <- Omega_tilde
  }
  
  fold.id <- sample(rep(seq(nfold),length = n))
  lambda.list <- array(list(),nfold);loglik <- array(list(),nfold)
  for(ii in 1:nfold){
    train.id <- fold.id != ii;valid.id <- fold.id == ii
    eval(parse(text = paste("x.train=x[", paste(rep(",", K), collapse = ""), "train.id]")))
    eval(parse(text = paste("x.valid=x[", paste(rep(",", K), collapse = ""), "valid.id]")))
    result <- Separate(x = x.train,val = x.valid,nlambda = nlambda[1],est.mode = est.mode,
                       rate = rate,normalize = normalize,lambda.list = NULL,
                       thr = thr,maxit = maxit)
    lambda.list[[ii]] <- result$lambda.list;loglik[[ii]] <- result$loglik
  }
  
  Lambda.list = Loglik <- array(list(),K)
  for(mode in 1:K){
    Lambda.list[[mode]] <- rowMeans(sapply(lambda.list,function(x) x[[mode]]))
    Loglik[[mode]] <- rowMeans(sapply(loglik, function(x) x[[mode]]))
  }
  
  lambda.list1 <- optimal.lambda.list(lambda = Lambda.list,loglik = Loglik,
                                      control.rate = control.rate,nlambda = nlambda[2])
  loglik1 <- array(list(),nfold)
  for(ii in 1:nfold){
    train.id <- fold.id != ii;valid.id <- fold.id == ii
    eval(parse(text = paste("x.train=x[", paste(rep(",", K), collapse = ""), "train.id]")))
    eval(parse(text = paste("x.valid=x[", paste(rep(",", K), collapse = ""), "valid.id]")))
    result <- Separate(x = x.train,val = x.valid,est.mode = est.mode,
                       rate = rate,normalize = normalize,lambda.list = lambda.list1,
                       thr = thr,maxit = maxit)
    loglik1[[ii]] <- result$loglik
  }
  
  lambda.best <- array(0,K)
  for(mode in 1:K){
    lambda.best[mode] <- lambda.list1[[mode]][which.max(rowMeans(sapply(loglik1,function(x) x[[mode]])))]
  }
  
  fit.result <- array(list(),K)
  for(mode_ind in 1:K){
    k = est.mode[mode_ind]
    Omega.list.sqrt = list()
    for (i in 1:K) {
      Omega.list.sqrt[[i]] = sqrtm(fit1[[i]])
    }
    # Calculate \tilde S_k
    S.array = array(0, c(m.vec[k], m.vec[k], n))
    Omega.list.sqrt[[k]] = diag(m.vec[k])
    for (i in 1:n) {
      d = 0
      eval(parse(text = paste("d=x[", paste(rep(",", K), collapse = ""), "i]")))
      Vi = k_unfold(as.tensor(ttl(as.tensor(d), Omega.list.sqrt,
                                  ms = 1:K
      )@data), m = k)@data
      S.array[, , i] = Vi %*% t(Vi)
    }
    S.mat = apply(S.array, c(1, 2), mean) * m.vec[k] / prod(m.vec) # \tilde S_k
    # fit model
    Out1 = glasso(S.mat, rho = lambda.best[mode_ind],
                  penalize.diagonal = FALSE, maxit = maxit, thr = thr)
    hat_Omega = as.matrix(Out1$wi)
    # normalization
    if (normalize) {
      hat_Omega = hat_Omega / norm(hat_Omega, type = "F")
    }
    fit.result[[mode_ind]] = hat_Omega
  }
  
  Result <- list(Omegahat = fit.result,
                 lambda.best = lambda.best,
                 loglik = loglik1,lambda.list = lambda.list1)
  return(Result)
}

Separate.penalty <- function(x,val = NULL,est.mode = NULL,Omega,
                             rate = 1e-3,nlambda = 20,lambda.list = NULL,
                             penalty = c("scad","atan","selo","mcp","exp","log","lq","sica","adapt"),
                             normalize = T,maxit = 1e3,thr = 1e-3, gamma = NULL){
  dimen = dim(x) # dimension of data set
  K = length(dimen) - 1 # order of tensor
  n = dimen[K + 1] # sample size of training set
  n.val <- dim(val)[K+1]
  nvars = prod(dimen) # number of variables
  m.vec = dimen[1:K] # dimension of each observation
  est.mode <- NULL
  if (is.null(est.mode) == TRUE) {
    est.mode = c(1:K)
  }
  if (is.null(gamma)){
    if (penalty == "scad") {
      gamma <- 3.7
    }else if (penalty == "mcp") {
      gamma <- 2
    }else if (penalty == "adapt") {
      gamma <- 1
    }else {
      gamma <- 0.01
    }
  }
  if(!is.null(lambda.list)){
    lambda.mat <- array(list(),K);loglik <- array(list(),K)
    for(mode in 1:K){
      Theta <- Omega[[mode]];lambda.mat[[mode]] <- array(list(),nlambda)
      for(i in 1:nlambda){
        lambda.mat[[mode]][[i]] <-  eval(parse(
          text =  paste0(
            penalty,
            "_deriv(Theta = Theta,lambda = lambda.list[[mode]][i], gamma = gamma)"
          )))
      }
    }
  }
  
  fit1 <- array(list(),K)
  for(mode in 1:K){
    if(n * prod(m.vec[-mode])<=((m.vec[mode]-1)/2*m.vec[mode])){
      Omega_tilde <- diag(m.vec[mode])
    }else{
      S.array <- array(0,c(m.vec[mode],m.vec[mode],n))
      for(i in 1:n){
        e <- 0
        eval(parse(text = paste("e=x[", paste(rep(",", K), collapse = ""), "i]"))) # assign the ith observation to e
        A <- k_unfold(as.tensor(e),mode)@data
        S.array[,,i] <- A %*% t(A)
      }
      S.mat <- apply(S.array,c(1,2),mean) / prod(m.vec[-mode]);Omega_tilde <- solve(S.mat)
      if (normalize) {
        Omega_tilde= Omega_tilde / norm(Omega_tilde, type = "F")}
    }
    fit1[[mode]] <- Omega_tilde
  }
  if(is.null(lambda.list)){
  lambda.list <- array(list(),K);lambda.mat <- array(list(),K);loglik <- array(list(),K)
  }
  if(!is.null(val)){
    Omega.list = list() # list of \tilde\Omega
    Omega.list.sqrt = list() # list of square root of \tilde\Omega
    for (k in 1:K) {
      Omega.list[[k]] = Omega[[k]]
      Omega.list.sqrt[[k]] = sqrtm(Omega.list[[k]])
    }
    Omega.sqrt.copy = Omega.list.sqrt
    
    for (mode_index in 1:K) {
      k = est.mode[mode_index]
      
      # Calculate \tilde S_k using the training set
      S.array = array(0, c(m.vec[k], m.vec[k], n))
      Omega.list.sqrt[[k]] = diag(m.vec[k]) # set \tilde\Omega_k to identity matrix
      for (i in 1:n){
        d = 0
        eval(parse(text = paste("d=x[", paste(rep(",", K), collapse = ""), "i]"))) # assign the ith observation to d
        Vi = k_unfold(as.tensor(ttl(as.tensor(d), Omega.list.sqrt,
                                    ms = 1:K
        )@data), m = k)@data
        S.array[, , i] = Vi %*% t(Vi)
      }
      S.mat = apply(S.array, c(1, 2), mean) * m.vec[k] / prod(m.vec) # \tilde S_k
      S.mat1 <- S.mat - diag(diag(S.mat))
      if(is.null(lambda.list[[mode_index]])){
      lambda.list[[mode_index]] <- exp(seq(log(rate * norm(S.mat1,"M")), log(norm(S.mat1,"M")),length.out = nlambda))
      Theta <- Omega[[mode_index]]
      for(i in 1:nlambda){
        lambda.mat[[mode_index]][[i]] <-  eval(parse(
          text =  paste0(
            penalty,
            "_deriv(Theta = Theta,lambda = lambda.list[[mode_index]][i], gamma = gamma)"
          )))
      }
      }
      
      # Calculate \tilde S_k using the validation set
      testS.array = array(0, c(m.vec[k], m.vec[k], n.val))
      for (i in 1:n.val) {
        d = 0
        eval(parse(text = paste("d=val[", paste(rep(",", K), collapse = ""), "i]")))
        Vi = k_unfold(as.tensor(ttl(as.tensor(d), Omega.list.sqrt,
                                    ms = 1:K
        )@data), m = k)@data
        testS.array[, , i] = Vi %*% t(Vi)
      }
      testS.mat = apply(testS.array, c(1, 2), mean) * m.vec[k] / prod(m.vec) # \tilde S_k
      Omega.list.sqrt[[k]] = Omega.sqrt.copy[[k]]
      
      # fit model with a sequence of lambdas
      lamk <- lambda.list[[mode_index]]
      lam.length = length(lamk)
      loglik2 = rep(0, lam.length)
      for (i in 1:lam.length){
        Out1 = glasso(S.mat, rho = lambda.mat[[mode_index]][[i]],
                      penalize.diagonal = FALSE, maxit = maxit, thr = thr)
        hat_Omega = Out1$wi;calc <- log(eigen(hat_Omega)$values)
        if(!is.complex(eigen(hat_Omega)$values)){
          loglik2[i] = -tr(testS.mat %*% hat_Omega) + sum(calc)
        }else{
          loglik2[i] <- NA
        }
      }
      loglik[[mode_index]] = loglik2
    }
  }
  
  lambda.best <- array(NA,K);mode.idx <- array(NA,K)
  for(mode in 1:K){
    mode.idx[mode] <- which.max(loglik[[mode]])
    lambda.best[mode] <- lambda.list[[mode]][which.max(loglik[[mode]])]
  }
  
  fit.result <- array(list(),K)
  for(mode_ind in 1:K){
    k = est.mode[mode_ind]
    Omega.list.sqrt = list()
    for (i in 1:K) {
      Omega.list.sqrt[[i]] = sqrtm(Omega[[i]])
    }
    # Calculate \tilde S_k
    S.array = array(0, c(m.vec[k], m.vec[k], n))
    Omega.list.sqrt[[k]] = diag(m.vec[k])
    for (i in 1:n) {
      d = 0
      eval(parse(text = paste("d=x[", paste(rep(",", K), collapse = ""), "i]")))
      Vi = k_unfold(as.tensor(ttl(as.tensor(d), Omega.list.sqrt,
                                  ms = 1:K
      )@data), m = k)@data
      S.array[, , i] = Vi %*% t(Vi)
    }
    S.mat = apply(S.array, c(1, 2), mean) * m.vec[k] / prod(m.vec) # \tilde S_k
    # fit model
    Out1 = glasso(S.mat, rho = lambda.mat[[mode_ind]][[mode.idx[mode_ind]]],
                  penalize.diagonal = FALSE, maxit = maxit, thr = thr)
    hat_Omega = as.matrix(Out1$wi)
    # normalization
    if (normalize) {
      hat_Omega = hat_Omega / norm(hat_Omega, type = "F")
    }
    fit.result[[mode_ind]] = hat_Omega
  }
  
  result <- list(Omegahat = fit.result,lambda.best = lambda.best,
                 loglik = loglik, lambda.list = lambda.list,
                 lambda.matrix = lambda.mat,penalty = penalty)
  return(result)
}

cv.Separate.penalty <- function(x,nfold = 5,est.mode = NULL, Omega,
                                rate = 1e-3,seed = 14,nlambda = c(20,10),
                                penalty = c("scad","atan","selo","mcp","exp","log","lq","sica","adapt"),
                                control.rate = 1e-2, normalize = T,
                                maxit = 1e3,thr = 1e-3, gamma = NULL){
  set.seed(seed)
  dimen = dim(x) # dimension of data set
  K = length(dimen) - 1 # order of tensor
  n = dimen[K + 1] # sample size of training set
  nvars = prod(dimen) # number of variables
  m.vec = dimen[1:K] # dimension of each observation
  
  if (is.null(est.mode) == TRUE) {
    est.mode = c(1:K)
  }
  
  fit1 <- array(list(),K)
  for(mode in 1:K){
    if(n * prod(m.vec[-mode])<=((m.vec[mode]-1)/2*m.vec[mode])){
      Omega_tilde <- diag(m.vec[mode])
    }else{
      S.array <- array(0,c(m.vec[mode],m.vec[mode],n))
      for(i in 1:n){
        e <- 0
        eval(parse(text = paste("e=x[", paste(rep(",", K), collapse = ""), "i]"))) # assign the ith observation to e
        A <- k_unfold(as.tensor(e),mode)@data
        S.array[,,i] <- A %*% t(A)
      }
      S.mat <- apply(S.array,c(1,2),mean) / prod(m.vec[-mode]);Omega_tilde <- solve(S.mat)
      if (normalize) {
        Omega_tilde= Omega_tilde / norm(Omega_tilde, type = "F")}
    }
    fit1[[mode]] <- Omega_tilde
  }
  if(is.null(gamma)){
    if (penalty == "scad") {
      gamma <- 3.7
    }else if (penalty == "mcp") {
      gamma <- 2
    }else if (penalty == "adapt") {
      gamma <- 1
    }else {
      gamma <- 0.01
    }
  }
  fold.id <- sample(rep(seq(nfold),length = n))
  lambda.list <- array(list(),nfold);loglik <- array(list(),nfold)
  for(ii in 1:nfold){
    train.id <- fold.id != ii;valid.id <- fold.id == ii
    eval(parse(text = paste("x.train=x[", paste(rep(",", K), collapse = ""), "train.id]")))
    eval(parse(text = paste("x.valid=x[", paste(rep(",", K), collapse = ""), "valid.id]")))
    result <- Separate.penalty(x = x.train,val = x.valid,
                               nlambda = nlambda[1],est.mode = est.mode,
                               rate = rate,normalize = normalize, penalty = penalty,
                               lambda.list = NULL,Omega = Omega,thr = thr,maxit = maxit)
    lambda.list[[ii]] <- result$lambda.list;loglik[[ii]] <- result$loglik
  }
  
  Lambda.list = Loglik <- array(list(),K)
  for(mode in 1:K){
    Lambda.list[[mode]] <- rowMeans(sapply(lambda.list,function(x) x[[mode]]))
    Loglik[[mode]] <- rowMeans(sapply(loglik, function(x) x[[mode]]))
  }
  
  lambda.list1 <- optimal.lambda.list(lambda = Lambda.list,loglik = Loglik,
                                      control.rate = control.rate,nlambda = nlambda[2])
  loglik1 <- array(list(),nfold)
  for(ii in 1:nfold){
    train.id <- fold.id != ii;valid.id <- fold.id == ii
    eval(parse(text = paste("x.train=x[", paste(rep(",", K), collapse = ""), "train.id]")))
    eval(parse(text = paste("x.valid=x[", paste(rep(",", K), collapse = ""), "valid.id]")))
    result <- Separate.penalty(x = x.train,val = x.valid,Omega = Omega,
                               est.mode = est.mode,penalty = penalty,
                               rate = rate,normalize = normalize,
                               lambda.list = lambda.list1,thr = thr,maxit = maxit)
    loglik1[[ii]] <- result$loglik
  }
  
  lambda.best <- array(0,K);mode.idx <- array(0,K)
  for(mode in 1:K){
    mode.idx[mode] <- which.max(rowMeans(sapply(loglik1,function(x) x[[mode]])))
    lambda.best[mode] <- lambda.list1[[mode]][which.max(rowMeans(sapply(loglik1,function(x) x[[mode]])))]
  }
  
  fit.result <- array(list(),K);lambda.best.mat <- array(list(),K)
  for(mode_ind in 1:K){
    k = est.mode[mode_ind]
    Omega.list.sqrt = list()
    for (i in 1:K) {
      Omega.list.sqrt[[i]] = sqrtm(Omega[[i]])
    }
    # Calculate \tilde S_k
    S.array = array(0, c(m.vec[k], m.vec[k], n))
    Omega.list.sqrt[[k]] = diag(m.vec[k])
    for (i in 1:n) {
      d = 0
      eval(parse(text = paste("d=x[", paste(rep(",", K), collapse = ""), "i]")))
      Vi = k_unfold(as.tensor(ttl(as.tensor(d), Omega.list.sqrt,
                                  ms = 1:K
      )@data), m = k)@data
      S.array[, , i] = Vi %*% t(Vi)
    }
    S.mat = apply(S.array, c(1, 2), mean) * m.vec[k] / prod(m.vec) # \tilde S_k
    # fit model
    lambda.best.mat[[mode_ind]] <-  eval(parse(
      text =  paste0(
        penalty,
        "_deriv(Theta = Omega[[mode_ind]],lambda = lambda.best[mode_ind], gamma = gamma)"
      )))
    Out1 = glasso(S.mat, rho = lambda.best.mat[[mode_ind]],
                  penalize.diagonal = FALSE, maxit = maxit, thr = thr)
    hat_Omega = as.matrix(Out1$wi)
    # normalization
    if (normalize) {
      hat_Omega = hat_Omega / norm(hat_Omega, type = "F")
    }
    fit.result[[mode_ind]] = hat_Omega
  }
  
  Result <- list(Omegahat = fit.result,
                 lambda.best = lambda.best,
                 loglik = loglik1,lambda.list = lambda.list1,
                 lambda.best.matrix = lambda.best.mat)
  return(Result)
}