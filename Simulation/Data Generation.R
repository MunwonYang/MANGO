library(MASS)
library(rTensor)
library(tensr)
library(ramify)
library(expm)
library(stats)
library(glasso)
library(Tlasso)

source("~/Basic Setting.R")

# Matrix Normal Data Generation
Generate.Normal <- function(n,dimen,Omega,seed){
  K <- length(Omega)
  Cov <- array(list(),K)
  stand <- array(list(),K)
  for(k in 1:K){
    Cov[[k]] <- solve(Omega[[k]])
  }
  
  alpha <- tr(Cov[[2]]) / dimen[2]
  
  beta <- sqrt(diag(Cov[[1]]))
  B <- diag(beta)
  
  Cov.Star <- array(list(),K)
  Cov.Star[[1]] <- solve(B) %*% Cov[[1]] %*% solve(B)
  Cov.Star[[2]] <- Cov[[2]] / alpha
  
  stand <- array(list(),K)
  for(k in 1:K){
    stand[[k]] <- sqrtm(Cov.Star[[k]])
  }
  
  set.seed(seed)
  X <- array(list(),n)
  Vax <- array(list(),n)
  
  x <- array(0,c(dimen,n))
  vax <- array(0,c(dimen,n))
  
  for(i in 1:n){
    Z <- array(rnorm(prod(dimen)),dimen)
    Z1 <- array(rnorm(prod(dimen)),dimen)
    X[[i]] <- atrans(Z,stand)
    Vax[[i]] <- atrans(Z1,stand)
    x[,,i] <- atrans(Z,stand)
    vax[,,i] <- atrans(Z1,stand)
  }
  
  Data <- list()
  Data$Separate.X <- X
  Data$Separate.Vax <- Vax
  Data$Cyclic.X <- x
  Data$Cyclic.Vax <- vax
  Omega.star <- array(list(),K)
  for(k in 1:K){
    Omega.star[[k]] <- solve(Cov.Star[[k]])
  }
  Data$Omega <- Omega.star
  return(Data)
}

# Matrix Data generated with Chi-squared distribution
Generate.Chi <- function(n,dimen,Omega,nu = NULL ,seed, i.i.d = TRUE,
                         variation = TRUE, nu.min = NULL, nu.max = NULL){
  K <- length(Omega)
  Cov <- array(list(),K)
  stand <- array(list(),K)
  for(k in 1:K){
    Cov[[k]] <- solve(Omega[[k]])
  }
  
  alpha <- tr(Cov[[2]]) / dimen[2]
  
  beta <- sqrt(diag(Cov[[1]]))
  B <- diag(beta)
  
  Cov.Star <- array(list(),K)
  Cov.Star[[1]] <- solve(B) %*% Cov[[1]] %*% solve(B)
  Cov.Star[[2]] <- Cov[[2]] / alpha
  
  stand <- array(list(),K)
  for(k in 1:K){
    stand[[k]] <- sqrtm(Cov.Star[[k]])
  }
  X <- array(list(),n); Vax <- array(list(),n)
  x <- array(0,c(dimen,n)); vax <- array(0,c(dimen,n))
  tail <- array(0,c(dimen[1],n))
  for(i in 1:n){
    set.seed(seed * i)
    if(variation == FALSE & !is.null(nu)){
      tail[,i] <- rep(nu/rchisq(1,nu),dimen[1])
    }else if(variation == TRUE & i.i.d == TRUE & !is.null(nu)){
      tail[,i] <- nu / rchisq(dimen[1],nu)
    }else if(variation == TRUE & i.i.d == FALSE & is.null(nu)){
      for(j in 1:dimen[1]){
        set.seed(j*i*seed)
        df <- runif(1,min = nu.min,max = nu.max)
        tail[j,i] <- sqrt(df/rchisq(1,df))
      }
    }
    Z <- array(rnorm(prod(dimen)),dimen)
    Z1 <- array(rnorm(prod(dimen)),dimen)
    X[[i]] <- diag(tail[,i])%*%atrans(Z,stand)
    Vax[[i]] <- diag(tail[,i])%*%atrans(Z1,stand)
    x[,,i] <- diag(tail[,i])%*%atrans(Z,stand)
    vax[,,i] <- diag(tail[,i])%*%atrans(Z1,stand)
  }
  
  Data <- list()
  Data$tail <- tail
  Data$Separate.X <- X
  Data$Separate.Vax <- Vax
  Data$Cyclic.X <- x
  Data$Cyclic.Vax <- vax
  Omega.star <- array(list(),K)
  for(k in 1:K){
    Omega.star[[k]] <- solve(Cov.Star[[k]])
  }
  Data$Omega <- Omega.star
  return(Data)
}

# Generating Contaminated Matrix data
Generate.Contaminated <- function(n,dimen,Omega,value,prob,seed,different = TRUE){
  if(sum(prob)!=1){
    stop("The Assigned probability is incorrect. The total sum of the probability should be 1")
  }
  if(sum(prob < 0) > 0){
    stop("The Assigned probability is incorrect. All probability should be same or larger than 0.")
  }
  if(length(value) != length(prob)){
    stop("Each random variable should have their assigned probability")
  }
  K <- length(Omega)
  Cov <- array(list(),K)
  stand <- array(list(),K)
  for(k in 1:K){
    Cov[[k]] <- solve(Omega[[k]])
  }
  
  alpha <- tr(Cov[[2]]) / dimen[2]
  
  beta <- sqrt(diag(Cov[[1]]))
  B <- diag(beta)
  
  Cov.Star <- array(list(),K)
  Cov.Star[[1]] <- solve(B) %*% Cov[[1]] %*% solve(B)
  Cov.Star[[2]] <- Cov[[2]] / alpha
  
  stand <- array(list(),K)
  for(k in 1:K){
    stand[[k]] <- sqrtm(Cov.Star[[k]])
  }
  X <- array(list(),n);Vax <- array(list(),n)
  x <- array(0,c(dimen,n)); vax <- array(0,c(dimen,n))
  tail <- array(0,c(dimen[1],n))
  for(i in 1:n){
    set.seed(i * seed)
    if(different == TRUE){
      tail[,i] <- sample(value,dimen[1],prob = prob,replace = TRUE)
    }else{
      tail[,i] <- rep(sample(value,1,prob = prob,replace = TRUE),dimen[1])
    }
    Z <- array(rnorm(prod(dimen)),dimen)
    Z1 <- array(rnorm(prod(dimen)),dimen)
    X[[i]] <- diag(tail[,i])%*%atrans(Z,stand)
    Vax[[i]] <- diag(tail[,i])%*%atrans(Z1,stand)
    x[,,i] <- diag(tail[,i])%*%atrans(Z,stand)
    vax[,,i] <- diag(tail[,i])%*%atrans(Z1,stand)
  }
  
  Data <- list()
  Data$tail <- tail
  Data$Separate.X <- X
  Data$Separate.Vax <- Vax
  Data$Cyclic.X <- x
  Data$Cyclic.Vax <- vax
  Omega.star <- array(list(),K)
  for(k in 1:K){
    Omega.star[[k]] <- solve(Cov.Star[[k]])
  }
  Data$Omega <- Omega.star
  return(Data)
}

# Generating Multivariate Laplace Distribution
Generate.Multivariate.Laplace <- function(n,dimen,Omega,seed = 14){
  K <- length(Omega)
  Cov <- array(list(),K)
  stand <- array(list(),K)
  for(k in 1:K){
    Cov[[k]] <- solve(Omega[[k]])
  }
  
  alpha <- tr(Cov[[2]]) / dimen[2]
  
  beta <- sqrt(diag(Cov[[1]]))
  B <- diag(beta)
  
  Cov.Star <- array(list(),K)
  Cov.Star[[1]] <- solve(B) %*% Cov[[1]] %*% solve(B)
  Cov.Star[[2]] <- Cov[[2]] / alpha
  
  stand <- array(list(),K)
  for(k in 1:K){
    stand[[k]] <- sqrtm(Cov.Star[[k]])
  }
  
  set.seed(seed)
  X <- array(list(),n);Vax <- array(list(),n)
  x <- array(0,c(dimen,n)); vax <- array(0,c(dimen,n))
  tail <- array(0,dimen[1]);tail <- abs(rexp(1) * rnorm(dimen[1]))
  
  for(i in 1:n){
    Z <- array(rnorm(prod(dimen)),dimen)
    Z1 <- array(rnorm(prod(dimen)),dimen)
    X[[i]] <- diag(tail)%*%atrans(Z,stand)
    Vax[[i]] <- diag(tail)%*%atrans(Z1,stand)
    x[,,i] <- diag(tail)%*%atrans(Z,stand)
    vax[,,i] <- diag(tail)%*%atrans(Z1,stand)
  }
  
  Data <- list()
  Data$tail <- tail
  Data$Separate.X <- X
  Data$Separate.Vax <- Vax
  Data$Cyclic.X <- x
  Data$Cyclic.Vax <- vax
  Omega.star <- array(list(),K)
  for(k in 1:K){
    Omega.star[[k]] <- solve(Cov.Star[[k]])
  }
  Data$Omega <- Omega.star
  return(Data)
}

Matrix.Identifiability <- function(Omega,thre = 1e-10){
  if(length(Omega)!=2){
    stop("The Data matrix is not Matrix Distribution.")
  }
  dimen <- c()
  for(i in 1:2){
    dimen <- c(dimen, length(diag(Omega[[i]])))
  }
  Cov <- array(list(),2)
  K <- length(Omega)
  for(k in 1:K){
    Cov[[k]] <- solve(Omega[[k]])
  }
  alpha <- tr(Cov[[2]]) / dimen[2]
  
  beta <- sqrt(diag(Cov[[1]]))
  B <- diag(beta)
  
  Cov.Tilde <- array(list(),2)
  Cov.Tilde[[1]] <- solve(B) %*% Cov[[1]] %*% solve(B)
  Cov.Tilde[[2]] <- Cov[[2]] / alpha
  
  Omega.Tilde <- array(list(),2)
  for(k in 1:2){
    Omega.Tilde[[k]] <- solve(Cov.Tilde[[k]])
    which <- abs(Omega.Tilde[[k]] - diag(diag(Omega.Tilde[[k]])))<=thre
    D <- diag(Omega.Tilde[[k]])
    Omega.Tilde[[k]][which] <- 0;Omega.Tilde[[k]] <- diag(D) + Omega.Tilde[[k]]
  }
  return(Omega.Tilde)
}