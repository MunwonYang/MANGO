library(Tlasso)
library(tensr)
library(expm)
library(rTensor)
library(glasso)

source("~/Heavy_tail/Basic Setting.R")

Separate.edit <- function(x,val,nlambda = c(50,20),
                          est.mode = NULL,rate = 1e-5, warm = TRUE, 
                          normalize = TRUE, thr = 1e-3, maxit = 1e3){
  dimen = dim(x) # dimension of dataset
  K = length(dimen) - 1 # order of tensor
  n = dimen[K + 1] # sample size of training set
  n_val = dim(val)[K + 1] # sample size of validation set
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
  
  lambda.list <- array(list(),K)
  loglik = loglik1 <- vector("list",K)
  for(mode in 1:K){
    if(norm(fit1[[mode]] - diag(diag(fit1[[mode]])),"M")==0){
      lambda.max <- 0.9
    }else{
      lambda.max <- 0.9 * norm(fit1[[mode]] - diag(diag(fit1[[mode]])),"M")
    }
    lambda.min <- rate * lambda.max
    lambda.list[[mode]] <-  exp(seq(log(lambda.min),log(lambda.max),length.out = nlambda[1]))
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
      
      # Calculate \tilde S_k using the validation set
      testS.array = array(0, c(m.vec[k], m.vec[k], n_val))
      for (i in 1:n_val) {
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
        loglik2[i] = -tr(testS.mat %*% hat_Omega) + sum(calc)
        if (loglik2[i] == Inf) {
          stop(paste("Infinite value! Please choose a smaller scale for mode", mode_index))
        }
        if (loglik2[i] == -Inf) {
          stop(paste("Negative infinite value! Please choose a larger scale for mode", mode_index))
        }
      }
      loglik[[mode_index]] = loglik2
    }
  }
  lambda.list1 <- optimal.lambda.list(lambda = lambda.list,
                                      loglik = loglik,
                                      nlambda = nlambda[2])
  lam.best = rep(0, K)
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
      
      # Calculate \tilde S_k using the validation set
      testS.array = array(0, c(m.vec[k], m.vec[k], n_val))
      for (i in 1:n_val) {
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
      lamk = lambda.list1[[mode_index]] # a sequence of candidates for lambda_k
      lam.length = length(lamk)
      loglik2 = rep(0, lam.length)
      for (i in 1:lam.length) {
        Out1 = glasso(S.mat, rho = lamk[i],
                      penalize.diagonal = FALSE, maxit = maxit, thr = thr)
        hat_Omega = Out1$wi;calc <- log(eigen(hat_Omega)$values)
        loglik2[i] = -tr(testS.mat %*% hat_Omega) + sum(calc)
        if (loglik2[i] == Inf) {
          stop(paste("Infinite value! Please choose a smaller scale for mode", mode_index))
        }
        if (loglik2[i] == -Inf) {
          stop(paste("Negative infinite value! Please choose a larger scale for mode", mode_index))
        }
      }
      ind = which.max(loglik2)
      lam.best[mode_index] = lamk[ind] # get the optimal lambda that maximizes the log-likelihood
      loglik1[[mode_index]] = loglik2
    }
    fit_result <- array(list(),K)
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
      Out1 = glasso(S.mat, rho = lam.best[mode_ind],
                    penalize.diagonal = FALSE, maxit = maxit, thr = thr)
      hat_Omega = as.matrix(Out1$wi)
      # normalization
      if (normalize) {
        hat_Omega = hat_Omega / norm(hat_Omega, type = "F")
      }
      fit_result[[mode_ind]] = hat_Omega
    }
  }
  
  result <- list(Omegahat = fit_result,
                 explore.lambda = lambda.list, explore.loglik = loglik,
                 final.lambda = lambda.list1, final.loglik = loglik1,
                 lambda = lam.best)
  return(result)
}

Separate.adapt <- function(x,val,Omega,est.mode = NULL,penalty = "adapt",
                           rate = 1e-5,nlambda = c(20,10),gam = 1,
                           normalize = TRUE, maxit = 1e3,
                           thr = 1e-3){
  dimen = dim(x) # dimension of dataset
  K = length(dimen) - 1 # order of tensor
  n = dimen[K + 1] # sample size of training set
  n_val = dim(val)[K + 1] # sample size of validation set
  nvars = prod(dimen) # number of variables
  m.vec = dimen[1:K] # dimension of each observation
  
  Lambda.list <- array(list(),K);weight.matrix <- array(list(),K)
  for(mode in 1:K){
    Theta <- Omega[[mode]]
    weight.matrix[[mode]] <-  eval(parse(
      text =  paste0(
        penalty,
        "_deriv(Theta = Theta, gam = gam)"
      )))
    if(norm(Omega[[mode]] - diag(diag(Omega[[mode]])),"M")==0){
      lambda.max <- 0.9
    }else{
      lambda.max <- 0.9 * norm(Omega[[mode]] - diag(diag(Omega[[mode]])),"M")
    }
    lambda.min <- rate * lambda.max
    Lambda.list[[mode]] <- exp(seq(log(lambda.min),log(lambda.max),length.out = nlambda[1]))
  }
  if (is.null(est.mode) == TRUE){
    est.mode = c(1:K)
  }
  loglik = loglik1 <- vector("list",K)
  if(!is.null(val)){
    Omega.list = list() # list of \tilde\Omega
    Omega.list.sqrt = list() # list of square root of \tilde\Omega
    for (k in 1:K) {
      Omega.list[[k]] = Omega[[k]]
      Omega.list.sqrt[[k]] = sqrtm(Omega.list[[k]])
    }
    Omega.sqrt.copy = Omega.list.sqrt
    for(mode_index in 1:K){
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
      
      # Calculate \tilde S_k using the validation set
      testS.array = array(0, c(m.vec[k], m.vec[k], n_val))
      for (i in 1:n_val) {
        d = 0
        eval(parse(text = paste("d=val[", paste(rep(",", K), collapse = ""), "i]")))
        Vi = k_unfold(as.tensor(ttl(as.tensor(d), Omega.list.sqrt,
                                    ms = 1:K
        )@data), m = k)@data
        testS.array[, , i] = Vi %*% t(Vi)
      }
      testS.mat = apply(testS.array, c(1, 2), mean) * m.vec[k] / prod(m.vec) # \tilde S_k
      Omega.list.sqrt[[k]] = Omega.sqrt.copy[[k]]
      
      lamk = Lambda.list[[mode_index]] # a sequence of candidates for lambda_k
      lam.length = length(lamk)
      loglik2 = rep(0, lam.length)
      
      for(i in 1:lam.length){
        Out1 <- glasso(S.mat,rho = lamk[i] * weight.matrix[[mode_index]],
                       penalize.diagonal = FALSE,
                       maxit = maxit, thr = thr)
        hat_Omega <- Out1$wi;calc <- log(eigen(hat_Omega)$values)
        loglik2[i] <- -tr(testS.mat %*% hat_Omega) + sum(calc)
        if (loglik2[i] == Inf) {
          stop(paste("Infinite value! Please choose a smaller scale for mode", mode_index))
        }
        if (loglik2[i] == -Inf) {
          stop(paste("Negative infinite value! Please choose a larger scale for mode", mode_index))
        }
      }
      loglik[[mode_index]] <- loglik2
    }
  }
  Lambda.list1 <- optimal.lambda.list(lambda = Lambda.list,
                                      loglik = loglik,
                                      nlambda = nlambda[2])
  lam.best <- rep(NA,K)
  if(!is.null(val)){
    Omega.list = list() # list of \tilde\Omega
    Omega.list.sqrt = list() # list of square root of \tilde\Omega
    for (k in 1:K) {
      Omega.list[[k]] = Omega[[k]]
      Omega.list.sqrt[[k]] = sqrtm(Omega.list[[k]])
    }
    Omega.sqrt.copy = Omega.list.sqrt
    for(mode_index in 1:K){
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
      
      # Calculate \tilde S_k using the validation set
      testS.array = array(0, c(m.vec[k], m.vec[k], n_val))
      for (i in 1:n_val) {
        d = 0
        eval(parse(text = paste("d=val[", paste(rep(",", K), collapse = ""), "i]")))
        Vi = k_unfold(as.tensor(ttl(as.tensor(d), Omega.list.sqrt,
                                    ms = 1:K
        )@data), m = k)@data
        testS.array[, , i] = Vi %*% t(Vi)
      }
      testS.mat = apply(testS.array, c(1, 2), mean) * m.vec[k] / prod(m.vec) # \tilde S_k
      Omega.list.sqrt[[k]] = Omega.sqrt.copy[[k]]
      
      lamk = Lambda.list1[[mode_index]] # a sequence of candidates for lambda_k
      lam.length = length(lamk)
      loglik2 = rep(0, lam.length)
      
      for(i in 1:lam.length){
        Out1 <- glasso(S.mat,rho = lamk[i] * weight.matrix[[mode_index]],
                       penalize.diagonal = FALSE,
                       maxit = maxit, thr = thr)
        hat_Omega <- Out1$wi;calc <- log(eigen(hat_Omega)$values)
        loglik2[i] <- -tr(testS.mat %*% hat_Omega) + sum(calc)
        if (loglik2[i] == Inf) {
          stop(paste("Infinite value! Please choose a smaller scale for mode", mode_index))
        }
        if (loglik2[i] == -Inf) {
          stop(paste("Negative infinite value! Please choose a larger scale for mode", mode_index))
        }
      }
      ind <- which.max(loglik2)
      lam.best[mode_index] <- lamk[ind]
      loglik1[[mode_index]] <- loglik2
    }
    fit_result <- array(list(),K)
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
      Out2 = glasso(S.mat, rho = lam.best[mode_ind] * weight.matrix[[mode_ind]],
                    penalize.diagonal = FALSE, maxit = maxit, thr = thr)
      hat_Omega = as.matrix(Out2$wi)
      # normalization
      if (normalize) {
        hat_Omega = hat_Omega / norm(hat_Omega, type = "F")
      }
      fit_result[[mode_ind]] = hat_Omega
    }
  }
  result <- list(Omegahat = fit_result,
                 lambda = lam.best,
                 explore.lambda = Lambda.list,
                 explore.loglik = loglik,
                 final.lambda = Lambda.list1,
                 final.loglik = loglik1,
                 weight = weight.matrix)
  return(result)
  
}