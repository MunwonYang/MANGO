library(MASS)
library(expm)
library(rTensor)
library(tensr)
library(glasso)
library(jointMeanCov)
library(stats)

source("~/Basic Setting.R")

Gemini.partial <- function(x,val = NULL,nlambda = c(50,20),
                           rate = 1e-5, warm = TRUE,
                           Omegatilde.list = NULL, scale.vec = NULL, 
                           normalize = TRUE, thr = 1.0e-3, maxit = 1e3){
  K <- length(dim(x))-1
  n <- dim(x)[K+1];n.val <- dim(val)[K+1]
  dimen <- dim(x)[1:K]
  for(i in 1:n.val){
    d=0
    eval(parse(text=paste('d=val[',paste(rep(',',K),collapse=''),'i]')))
    if(sum(dim(d)!=dimen)!=0){
      stop("Dimension of Validation data is different with Training data.")
    }
  }
  fit.x = fit.vax <- array(list(),K)
  for(mode in 1:K){
    order <- seq(1:length(dim(x)));order[mode] <- NA
    order <- as.numeric(na.omit(c(mode,order)))
    corr <- aperm(x,order);corr.vax <- aperm(val,order)
    fit.x[[mode]] <- corr.est(corr);fit.vax[[mode]] <- corr.est(corr.vax)
  }
  lambda.list <- array(list(),K)
  for(mode in 1:K){
    if(norm(fit.x[[mode]] - diag(diag(fit.x[[mode]])),"M") == 0){
      lambda.max <- 0.9
    }else{
      lambda.max <- 0.9 * norm(fit.x[[mode]] - diag(diag(fit.x[[mode]])),"M")
    }
    lambda.min <- rate * lambda.max
    lambda.list[[mode]] <-  exp(seq(log(lambda.min),log(lambda.max),length.out = nlambda[1]))
  }
  loglik <- vector("list",K)
  for(mode in 1:K){
    lambda <- lambda.list[[mode]];lam.length <- length(lambda)
    log.lik <- array(0,lam.length)
    for(i in 1:lam.length){
      lam <- lambda[i]
      out <- glasso(fit.x[[mode]],rho = lam, 
                    penalize.diagonal = FALSE, 
                    maxit = maxit, thr = thr)$wi
      calc <- log(eigen(out)$values)
      log.lik[i] <- tr(out%*%fit.vax[[mode]]) - sum(calc)
    }
    loglik[[mode]] <- log.lik
  }
  lambda.list1 <- optimal.lambda.list(lambda = lambda.list,type = "convex",
                                      loglik = loglik, nlambda = nlambda[2])
  loglik1 <- vector("list",K)
  lambda.best <- array(NA,K)
  for(mode in 1:K){
    lambda <- lambda.list1[[mode]];lam.length <- length(lambda)
    log.lik <- array(0,lam.length)
    for(i in 1:lam.length){
      lam <- lambda[i]
      out <- glasso(fit.x[[mode]],rho = lam, 
                    penalize.diagonal = FALSE, 
                    maxit = maxit, thr = thr)$wi
      calc <- log(eigen(out)$values)
      log.lik[i] <- tr(out%*%fit.vax[[mode]]) - sum(calc)
    }
    lambda.best[mode] <- lambda[which.min(log.lik)]
    loglik1[[mode]] <- log.lik
  }
  fit.output <- array(list(),K)
  fit.result <- array(list(),K)
  fit.sd.row <- array(list(),K)
  fit.inv.sd.row <- array(list(),K)
  fit.dimen.setting <- array(list(),K)
  for(mode in 1:K){
    order <- seq(1:length(dim(x)));order[mode] <- NA
    order <- as.numeric(na.omit(c(mode,order)))
    output <- Gemini.edit(aperm(x,order),lambda = lambda.best[mode])
    fit.output[[mode]] <- output$Pre
    fit.result[[mode]] <- output$Corr
    fit.sd.row[[mode]] <- output$sd.row
    fit.inv.sd.row[[mode]] <- output$inv.sd.row
    fit.dimen.setting[[mode]] <- dim(aperm(x,order))[1:K]
    if(normalize){
      fit.output[[mode]] <- fit.output[[mode]] / norm(fit.output[[mode]],"F")
    }
  }
  result <- list()
  result$Omegahat <- fit.output;result$lambda <- lambda.best
  result$corr = fit.result;result$inv.sd.row <- fit.inv.sd.row
  result$explore.lambda <- lambda.list;result$final.lambda <- lambda.list1
  result$explore.loglik <- loglik;result$final.loglik <- loglik1
  return(result)
}

Gemini.adapt <- function(x,val = NULL, input,inv.sd.row,
                         rate = 1e-5,bias = 1,
                         nlambda = c(20,10),normalize = TRUE,
                         gam = 1, thr = 1e-3, maxit = 1e3,
                         scale.vec = NULL){
  K <- length(dim(x))-1
  n <- dim(x)[K+1];n.val <- dim(val)[K+1]
  dimen <- dim(x)[1:K]
  for(i in 1:n.val){
    d=0
    eval(parse(text=paste('d=val[',paste(rep(',',K),collapse=''),'i]')))
    if(sum(dim(d)!=dimen)!=0){
      stop("Dimension of Validation data is different with Training data.")
    }
  }
  fit.x = fit.vax <- array(list(),K)
  for(mode in 1:K){
    order <- seq(1:length(dim(x)));order[mode] <- NA
    order <- as.numeric(na.omit(c(mode,order)))
    corr <- aperm(x,order);corr.vax <- aperm(val,order)
    fit.x[[mode]] <- corr.est(corr);fit.vax[[mode]] <- corr.est(corr.vax)
  }
  lambda.list <- array(list(),K)
  weight.matrix <- array(list(),K)
  for(mode in 1:K){
    weight.matrix[[mode]] <- adapt_deriv(input[[mode]],gam = gam) * bias
    if(norm(input[[mode]] - diag(diag(input[[mode]])),"M") == 0){
      lambda.max <- 0.9
    }else{
      lambda.max <- 0.9 * norm(input[[mode]] - diag(diag(input[[mode]])),"M");
    }
    lambda.min <- rate*lambda.max
    lambda.list[[mode]] <- exp(seq(log(lambda.min),log(lambda.max),length.out = nlambda[1]))
  }
  loglik <- vector("list",K)
  for(mode in 1:K){
    lambda <- lambda.list[[mode]];lam.length <- length(lambda)
    log.lik <- array(0,lam.length)
    for(i in 1:lam.length){
      lam <- lambda[i]
      out <- glasso(input[[mode]],rho = lam * weight.matrix[[mode]], 
                    penalize.diagonal = FALSE, 
                    maxit = maxit, thr = thr)$wi
      calc <- log(eigen(out)$values)
      log.lik[i] <- tr(out%*%fit.vax[[mode]]) - sum(calc)
    }
    loglik[[mode]] <- log.lik
  }
  lambda.list1 <- optimal.lambda.list(lambda = lambda.list,type = "convex",
                                      loglik = loglik,nlambda = nlambda[2],
                                      adjust = 10)
  loglik1 <- vector("list",K)
  lambda.best <- array(NA,K)
  for(mode in 1:K){
    lambda <- lambda.list1[[mode]];lam.length <- length(lambda)
    log.lik <- array(0,lam.length)
    for(i in 1:lam.length){
      lam <- lambda[i]
      out <- glasso(fit.x[[mode]],rho = lam * weight.matrix[[mode]], 
                    penalize.diagonal = FALSE, 
                    maxit = maxit, thr = thr)$wi
      calc <- log(eigen(out)$values)
      log.lik[i] <- tr(out%*%fit.vax[[mode]]) - sum(calc)
    }
    lambda.best[mode] <- lambda[which.min(log.lik)]
    loglik1[[mode]] <- log.lik
  }
  
  dimen.edit <- array(list(),K)
  fit.result <- array(list(),K)
  for(mode in 1:K){
    order <- seq(1:length(dim(x)));order[mode] <- NA
    order <- as.numeric(na.omit(c(mode,order)))
    dimen.edit[[mode]] <- dim(aperm(x,order))[1:K]
    out1 <- glasso(input[[mode]],rho = lambda.best[mode] * weight.matrix[[mode]], 
                   penalize.diagonal = FALSE, 
                   maxit = maxit, thr = thr)$wi
    fit.result[[mode]] <- dimen.edit[[mode]][2] * t(t(out1 * inv.sd.row[[mode]]) * inv.sd.row[[mode]])
  }
  result <- list()
  result$weight = weight.matrix;result$Omegahat <- fit.result
  result$lambda.best <- lambda.best;result$explore.lambda <- lambda.list
  result$explore.loglik <- loglik;result$final.lambda <- lambda.list1;result$final.loglik <- loglik1
  return(result)
}