# Covariance Setting
library(Tlasso)
library(MASS)
library(expm)
library(rTensor)
library(tensr)
library(glasso)

# Triangle Covariance
TR <- function(p){
  S <- ChainOmega(p, sd = p*  100, norm.type = 2)
  return(S)
}

# Auto regressive precision
AR <- function(ratio,p){
  S <- diag(p)
  for(i in 1:p){
    for (j in 1:p){
      S[i,j] = ratio^abs(i-j)
    }
  }
  return(S)
}

# Compound symmetry precision
CS <- function(p){
  S <- diag(p)
  for (i in 1:p){
    for (j in 1:p){
      if(i != j){
        S[i,j] = 0.6
      }
    }
  }
  return(S)
}

unfold <- function(X,mode){
  dimen = dim(X)
  order = seq(1,length(dimen))
  order[mode] = NA
  order <- as.numeric(na.omit(c(mode,order)))
  A <- array(aperm(X,order),c(dimen[mode],prod(dimen)/dimen[mode]))
  reutrn(A)
}

# Lambda selecting rule
optimal.lambda.list <- function(lambda,loglik,type = "concave",
                                nlambda = 20, control.rate =1e-4,
                                adjust = 3){
  K <- length(lambda)
  lambda.idx <- array(list(),K)
  lambda.list <- vector("list",K)
  if(type == "concave"){
    for(mode in 1:K){
      lambda.idx[[mode]] <- which.max(loglik[[mode]])
      if(lambda.idx[[mode]] == 1 | lambda.idx[[mode]] == 2){
        lambda.min <- control.rate * lambda[[mode]][[lambda.idx[[mode]]+adjust]]
        lam <- seq(log(lambda.min),log(lambda[[mode]][[lambda.idx[[mode]]+adjust]]),length.out = nlambda)
        lambda.list[[mode]] <- exp(lam)
      }else if(lambda.idx[[mode]] == length(lambda[[mode]]) | lambda.idx[[mode]] == (length(lambda[[mode]])-1)){
        lambda.max <- lambda[[mode]][[lambda.idx[[mode]]-adjust]]/control.rate
        lam <- seq(log(lambda[[mode]][[lambda.idx[[mode]]-adjust]]),log(lambda.max),length.out = nlambda)
        lambda.list[[mode]] <- exp(lam)
      }else{
        lambda.min <- lambda[[mode]][lambda.idx[[mode]]-1]
        lambda.max <- lambda[[mode]][lambda.idx[[mode]]+1]
        lambda.list[[mode]] <- seq(lambda.min,lambda.max,length.out = nlambda)
      }
    }
  }else{
    for(mode in 1:K){
      lambda.idx[[mode]] <- which.min(loglik[[mode]])
      if(lambda.idx[[mode]] == 1 | lambda.idx[[mode]] == 2){
        lambda.min <- control.rate * lambda[[mode]][[lambda.idx[[mode]]+adjust]]
        lam <- seq(log(lambda.min),log(lambda[[mode]][[lambda.idx[[mode]]+adjust]]),length.out = nlambda)
        lambda.list[[mode]] <- exp(lam)
      }else if(lambda.idx[[mode]] == length(lambda[[mode]]) | lambda.idx[[mode]] == (length(lambda[[mode]])-1)){
        lambda.max <- lambda[[mode]][[lambda.idx[[mode]]-adjust]]/control.rate
        lam <- seq(log(lambda[[mode]][[lambda.idx[[mode]]-adjust]]),log(lambda.max),length.out = nlambda)
        lambda.list[[mode]] <- exp(lam)
      }else{
        lambda.min <- lambda[[mode]][lambda.idx[[mode]]-1]
        lambda.max <- lambda[[mode]][lambda.idx[[mode]]+1]
        lambda.list[[mode]] <- seq(lambda.min,lambda.max,length.out = nlambda)
      }
    }
  }
  return(lambda.list)
}

adapt_deriv_weight <- function(Theta,gamma = 1){
  Theta <- (abs(Theta) + 1e-2)
  weight<- Theta^(-gamma)
  return(weight)
}


# Estimated the moment based correlation 
corr.est <- function(X){
  if(length(dim(X))==2){
    X <- array(X,c(dim(X),1))
  }
  XX <- 0;X2 <- 0
  dimen <- dim(X);n <- dimen[3]
  for(i in 1:n){
    XX <- XX + X[,,i]%*%t(X[,,i])
    X2 <- X2 + rowSums(X[,,i]^2)
  }
  X2sqrt <- sqrt(X2)
  Gamma.hat.B <- t(t(XX/X2sqrt)/X2sqrt)
  return(Gamma.hat.B)
}

Matrix.Surrogate <- function(x, mode = 1 ,type = c("Constant","Matrix","Original")){
  n <- dim(x)[length(dim(x))];K <- length(dim(x)) - 1;dimen <- dim(x)[1:K]
  if(type == "Constant"){
    gamma <- array(0,n)
    for(i in 1:n){
      d <- 0
      eval(parse(text = paste0("d=x[", paste0(rep(",", K), collapse = ""), "i]")))
      gamma[i] <- norm(d,"F") / sqrt(prod(dimen))
    }
  }else if(type == "Matrix"){
    gamma <- array(0,c(dimen[mode],n))
    for(j in 1:dimen[mode]){
      for(i in 1:n){
        d <- 0
        eval(parse(text = paste0("d = x[",paste0(rep(",",mode-1),"j",paste0(rep(",",K -(mode - 1)),collapse = ""),"i]"))))
        gamma[j,i] <- norm(d,"2")/sqrt(prod(dimen[setdiff(1:K,mode)]))
      }
    }
  }else{
    gamma <- array(1,n)
  }
  return(gamma)
}

Matrix.tail.test <- function(x,type = c("type1","type2"),
                             alternative = c("two.sided", "less", "greater")){
  K <- length(dim(x)) - 1
  if(K !=2){
    stop("The data needs to be a matrix.")
  }
  dimen <- dim(x)[1:K];n <- dim(x)[K+1]
  mat <- Matrix.Surrogate(x,type = "Matrix");con <- Matrix.Surrogate(x,type = "Constant")
  if(length(type)!=1){
    type <- "type1"
  }
  if(type == "type1"){
    mat <- Matrix.Surrogate(x,type = "Matrix")
    a <- mean(mat^2);b <- mean(mat^4)
    z.mat <- (mat^2 - a) / sqrt(2 * b / (dimen[2]+2))
    test.stat <- sum((z.mat - mean(z.mat))^2)
    p.val <- pchisq(test.stat,df = dimen[1]*n-1,lower.tail = F)
    alternative <- "greater"
  }else{
    est <- mean(apply(sweep(mat^2,2,con^2)^2,2,mean))
    p1 <- dimen[1];p2 <- dimen[2];p <- prod(dimen)
    E1 <- mean(con^4) * (p2/ (p+2));Var.con <- sum((con^4 - mean(con^4))^2) /(n-1)
    C1 <- (1 + 2/p) * (1 + 4/p) * (1 + 6/p);C2 <- ((p+2)/p)^2;C <- C1 - C2
    a1 <- (1 / p1)^2;a2 <- Var.con - (p/(p+2)) * C * mean(con^4)
    Var1 <- a1/C1 * a2
    E2 <- (2*dimen[1] - 2) / dimen[2]
    D1 <- (dimen[1]-1)^2/dimen[1]
    D2 <- (12 * (dimen[2]+4)) / dimen[2]^3
    D3 <- (dimen[1]-3) / (dimen[1]-1);D4 <- 4 / dimen[2]^2
    Var2 <- D1 * (D2 - D3 * D4)
    x.bar <- E1 * E2
    var.est <- E1^2 * Var2 + Var1 * E2^2 + Var1 * Var2
    test.stat <- (est - x.bar) / (sqrt(var.est/n))
    if(length(alternative)!=1){
      alternative <- "two.sided"
      if(pnorm(test.stat)<= 0.5){
        p.val <- pnorm(test.stat) * 2
      }else{
        p.val <- pnorm(test.stat,lower.tail = F) * 2
      }
    }else{
      if(alternative == "two.sided"){
        if(pnorm(test.stat)<= 0.5){
          p.val <- pnorm(test.stat) * 2
        }else{
          p.val <- pnorm(test.stat,lower.tail = F) * 2
        }
      }else if(alternative == "greater"){
        p.val <- pnorm(test.stat,lower.tail = F)
      }else{
        p.val <- pnorm(test.stat)
      }
    }
  }
  if(type == "type1"){
    est <- z.mat;estimate  <- NULL
  }else{
    if(is.null(x.bar)){
      estimate = var.est
    }else{
      estimate <-  data.frame(mean = x.bar,variance = var.est)
    }
  }
  value <- list(statistic = test.stat,
                p.value = p.val, sample.estimate = est,
                population.estimate = estimate,
                alternative = alternative)
  return(value)
}

Matrix.power <- function(x,type = c("type1","type2"),alpha = 0.05,
                         alternative = c("one.sided","two.sided")){
  K <- length(dim(x)) - 1
  if(K !=2){
    stop("The data needs to be a matrix.")
  }
  dimen <- dim(x)[1:K];n <- dim(x)[K+1]
  mat <- Matrix.Surrogate(x,type = "Matrix");con <- Matrix.Surrogate(x,type = "Constant")
  if(length(type)!=1){
    type <- "type1"
  }
  if(type == "type1"){
    mat <- Matrix.Surrogate(x,type = "Matrix")
    a <- mean(mat^2);b <- mean(mat^4)
    z.mat <- (mat^2 - a) / sqrt(2 * b / (dimen[2]+2))
    test.stat <- sum((z.mat - mean(z.mat))^2);power <- pchisq(test.stat,df = dimen[1]*n-1)
  }else{
    est <- mean(apply(sweep(mat^2,2,con^2)^2,2,mean))
    p1 <- dimen[1];p2 <- dimen[2];p <- prod(dimen)
    E1 <- mean(con^4) * (p2/ (p+2));Var.con <- sum((con^4 - mean(con^4))^2) /(n-1)
    C1 <- (1 + 2/p) * (1 + 4/p) * (1 + 6/p);C2 <- ((p+2)/p)^2;C <- C1 - C2
    a1 <- (1 / p1)^2;a2 <- Var.con - (p/(p+2)) * C * mean(con^4)
    Var1 <- a1/C1 * a2
    E2 <- (2*dimen[1] - 2) / dimen[2]
    D1 <- (dimen[1]-1)^2/dimen[1]
    D2 <- (12 * (dimen[2]+4)) / dimen[2]^3
    D3 <- (dimen[1]-3) / (dimen[1]-1);D4 <- 4 / dimen[2]^2
    Var2 <- D1 * (D2 - D3 * D4)
    x.bar <- E1 * E2
    var.est <- E1^2 * Var2 + Var1 * E2^2 + Var1 * Var2
    test.stat <- (est - x.bar) / (sqrt(var.est/n))
    if(length(alternative)!=1){
      alternative <- "two.sided"
      power <- pnorm(qnorm(alpha/2) + abs(test.stat))
    }else{
      if(alternative == "one.sided"){
        if(est>=x.bar){
          power <- pnorm(qnorm(alpha) + test.stat)
        }else{
          power <- 1 - pnorm(qnorm(1-alpha) + test.stat)
        }
      }else{
        power <- pnorm(qnorm(alpha/2) + abs(test.stat))
      }
    }
  }
  if(type == "type1"){
    est <- z.mat;estimate <- NULL
  }else{
    if(is.null(x.bar)){
      estimate = var.est
    }else{
      estimate <-  data.frame(mean = x.bar,variance = var.est)
    }
  }
  value <- list(statistic = test.stat,
                power = power, sample.estimate = est,
                population.estimate = estimate,
                alternative = alternative)
  return(value)
}

Hypothesis.test.parameter <- function(n,dimen,N = 1000,
                                      seed.value = 124, ncores = 96){
  Omega <- list(diag(dimen[1]),diag(dimen[2]))
  T.star <- foreach(iter = 1:N,
                    .combine = rbind, .packages = c("MASS","rTensor","tensr","expm","stats"))%dopar%{
                      source("~/Basic Setting.R")
                      source("~/Hypothesis Test/Basic Setting for Hypothesis test.R")
                      seed <- iter * seed.value 
                      X <- Generate.Normal(n = n,dimen = dimen,Omega = Omega,seed = seed)
                      mat <- Matrix.Surrogate(X$x,type = "Matrix");con <- Matrix.Surrogate(X$x,type = "Constant")
                      a <- mean((sweep(mat^2,2,con^2,"/")-1)^2)
                      return(a)
                    }
  result <- list(mu = mean(T.star), sigma = sd(T.star))
  return(result)
}


# Find the index for jth element in tensor structure
indexing <- function(a, dim){
  K = length(dim)
  s <- rep(1,K+1)
  j <- rep(1,K)
  s[K+1] <- a
  for(i in seq(K,2,-1)){
    s[i] <- s[i+1] %%prod(dim[1:(i-1)])
  }
  j[1] <- s[2]
  if(j[1] == 0){
    j[1] <- dim[1]
  }
  for(i in 2:K){
    j[i] <- ceiling(s[i+1] / prod(dim[1:(i-1)]))
    if(j[i] == 0){
      j[i] <- dim[i]
    }
  }
  j
}


# Express Group of index
index<- function(ix,dim){
  n <- length(ix);K <- length(dim)
  index.set <- array(0,c(n,K))
  for(i in 1:n){
    index.set[i,] <- indexing(ix[i],dim = dim)
  }
  return(index.set)
}

# Finding the vectorize index from tensor index when the tensor data vectorize
index.vectorize <- function(idx,dim){
  if(length(idx)!=length(dim)){
    stop("The length of index and dimension list needs to be identical.")
  }
  K <- length(dim)
  A <- c()
  for(i in 1:(K-1)){
    A <- c(A,prod(dim[1:i]))
  }
  loc <- 0
  for(j in 2:K){
    loc <- loc + (idx[j]-1) * A[(j-1)]
  }
  loc <- loc + idx[1]
  return(loc)
}

# Express the group of index when tensor data vectorize
index.vec <- function(ix,dim){
  if(ncol(ix)!=length(dim)){
    stop("The length of index and dimension list needs to be identical.")
  }
  n <- nrow(ix)
  index <- c()
  for(i in 1:n){
    index <- c(index,index.vectorize(ix[i,],dim = dim))
  }
  return(index)
}