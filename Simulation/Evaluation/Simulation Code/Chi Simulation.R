library(MASS)
library(expm)
library(rTensor)
library(tensr)
library(ramify)
library(stats)
library(glasso)
library(Tlasso)
library(doParallel)
library(foreach)
library(stats)
library(jointMeanCov)

source("~/Basic Setting.R")
source("~/Simulation Evaluation.R")
source("~/Simulation/Data Generation.R")
source("~/Simulation/Evaluation/Simulation Code/Methods/Separate.R")
source("~/Simulation/Evaluation/Simulation Code/Methods/Gemini.R")
source("~/Simulation/Evaluation/Simulation Code/Methods/Tlasso.R")

Separate.simul <- function(n, dimen, gam = 1,
                           penalty = 'adapt',Omega,scale = NULL,
                           seed = 14, norm.off.diag = FALSE,
                           rate = c(1e-5,1e-6),
                           Heavy = "Matrix",nu = NULL,i.i.d = TRUE,
                           variation = TRUE, nu.min = NULL, nu.max = NULL,
                           Run = 100, sparse = c(1,2)){
  Separate <- foreach(j = 1:Run, .combine =rbind,
                      .packages = c("tensr","ramify","stats",'rTensor','expm',
                                    "Tlasso","glasso","MASS"))%dopar%{
                                      source("~/Basic Setting.R")
                                      source("~/Simulation Evaluation.R")
                                      source("~/Simulation/Evaluation/Simulation Code/Methods/Separate.R")
                                      source("~/Simulation/Data Generation.R")
                                      Data <- Generate.Chi(n = n, dimen = dimen, Omega = Omega,
                                                           nu = nu, seed = seed*j, i.i.d = i.i.d,
                                                           variation = variation, nu.min = nu.min,
                                                           nu.max = nu.max)
                                      X <- Data$Cyclic.X;Vax <- Data$Cyclic.Vax
                                      K <- length(dimen);n.train <- dim(X)[K+1];n.val <- dim(Vax)[K+1]
                                      if(Heavy == "Matrix"){
                                        gamma <- array(0,c(dimen[1],n.train));gamma.vax <- array(0,c(dimen[1],n.val))
                                        for(j in 1:dimen[1]){
                                          for(i in 1:n.train){
                                            gamma[j,i] <- norm(X[j,,i],"2")/sqrt(dimen[2])
                                          }
                                          for(i in 1:n.val){
                                            gamma.vax[j,i] <- norm(Vax[j,,i],"2")/sqrt(dimen[2])
                                          }
                                        }
                                        x <- array(0,c(dimen,n.train));val <- array(0,c(dimen,n.val))
                                        for(i in 1:n.train){
                                          x[,,i] <- solve(diag(gamma[,i])) %*% X[,,i]
                                        }
                                        for(i in 1:n.val){
                                          val[,,i] <- solve(diag(gamma.vax[,i])) %*% Vax[,,i]
                                        }
                                      }else if(Heavy == "Constant"){
                                        gamma <- array(0,n.train);gamma.vax <- array(0,n.val)
                                        x <- array(0,c(dimen,n.train));val <- array(0,c(dimen,n.val))
                                        for(i in 1:n.train){
                                          gamma[i] <- norm(X[,,i],"F")/sqrt(prod(dimen))
                                          x[,,i] <- X[,,i] / gamma[i]
                                        }
                                        for(i in 1:n.val){
                                          gamma.vax[i] <- norm(Vax[,,i],"F")/sqrt(prod(dimen))
                                          val[,,i] <- Vax[,,i] / gamma.vax[i]
                                        }
                                      }else{
                                        x <- X;val <- Vax
                                      }
                                      fit <- Separate.edit(x = x, val = val, rate = rate[1])
                                      if(penalty == 'adapt'){
                                        fit.adapt <- Separate.adapt(x = x, val = val, rate = rate[2],
                                                                    Omega = fit$Omegahat)
                                        out <- Matrix.simulation.summary(Omega.hat.list = fit.adapt$Omegahat,
                                                                         Omega.true.list = Omega,
                                                                         norm.off.diag = norm.off.diag)
                                      }else{
                                        out <- Matrix.simulation.summary(Omega.hat.list = fit$Omegahat,
                                                                         Omega.true.list = Omega,
                                                                         norm.off.diag = norm.off.diag)
                                      }
                                      result <- list(Avg.F = out$av.error.f*100, Avg.Max = out$av.error.max*100,
                                                     Error.F = out$error.f*100, Error.Max = out$error.max*100,
                                                     Avg.Tpr = out$av.tpr*100, Avg.Fpr = out$av.fpr*100,
                                                     Avg.Tnr = out$av.tnr*100,
                                                     Tpr = out$tpr*100, Fpr = out$fpr*100,
                                                     Tnr = out$tnr * 100,
                                                     Optimality = fit$optimal)
                                      result
                                    }
  Result <- as.data.frame(Separate)
  Error.F <- matrix(unlist(Result$Error.F),Run,byrow = TRUE)
  Error.Max <- matrix(unlist(Result$Error.Max),Run,byrow = TRUE)
  TPR <- matrix(unlist(Result$Tpr),Run, byrow = TRUE)
  FPR <- matrix(unlist(Result$Fpr),Run, byrow = TRUE)
  K = length(dimen)
  if(sum(sparse == seq(1:K))==K){
    Outcome <- list(Avg.Error.F.mean = round(mean(unlist(Result$Avg.F)),2), Avg.Error.F.sd =   round(sd(unlist(Result$Avg.F))/sqrt(Run),2),
                    Avg.Error.Max.mean = round(mean(unlist(Result$Avg.Max)),2), Avg.Error.Max.sd =   round(sd(unlist(Result$Avg.Max))/sqrt(Run),2),
                    Avg.TPR.mean = round(mean(unlist(Result$Avg.Tpr)),2), Avg.TPR.sd =   round(sd(unlist(Result$Avg.Tpr))/sqrt(Run),2),
                    Avg.FPR.mean = round(mean(unlist(Result$Avg.Fpr)),2), Avg.FPR.sd =   round(sd(unlist(Result$Avg.Fpr))/sqrt(Run),2),
                    Error.F.mean = round(colMeans(Error.F),2), Error.F.sd = round(apply(Error.F,2,sd)/sqrt(Run),2),
                    Error.Max.mean = round(colMeans(Error.Max),2), Error.Max.sd = round(apply(Error.Max,2,sd)/sqrt(Run),2),
                    TPR.mean = round(colMeans(TPR),2), TPR.sd = round(apply(TPR,2,sd)/sqrt(Run),2),
                    FPR.mean = round(colMeans(FPR),2), FPR.sd = round(apply(FPR,2,sd)/sqrt(Run),2))
  }
  if(sum(sparse == seq(1:K))!=K){
    Outcome <- list(Error.F.mean = round(colMeans(Error.F),2)[sparse], Error.F.sd = round(apply(Error.F,2,sd)/sqrt(Run),2)[sparse],
                    Error.Max.mean = round(colMeans(Error.Max),2)[sparse], Error.Max.sd = round(apply(Error.Max,2,sd)/sqrt(Run),2)[sparse],
                    TPR.mean = round(colMeans(TPR),2)[sparse], TPR.sd = round(apply(TPR,2,sd)/sqrt(Run),2)[sparse],
                    FPR.mean = round(colMeans(FPR),2)[sparse], FPR.sd = round(apply(FPR,2,sd)/sqrt(Run),2)[sparse])
  }
  return(Outcome)
}

Cyclic.simul <- function(n, dimen, Omega,scale = NULL,
                         seed = 14, norm.off.diag = FALSE,
                         gam = 1, rate = c(1e-5,1e-6),
                         penalty = 'adapt',Heavy = "Matrix",
                         nu = NULL,i.i.d = TRUE,thres = 1e-5,t = 5,
                         variation = TRUE, nu.min = NULL, nu.max = NULL,
                         Run = 100, sparse = c(1,2)){
  Cyclic <- foreach(j = 1:Run, .combine =rbind,
                    .packages = c("tensr","ramify","stats",'rTensor','expm',
                                  "Tlasso","glasso","MASS"))%dopar%{
                                    source("~/Basic Setting.R")
                                    source("~/Simulation Evaluation.R")
                                    source("~/Simulation/Evaluation/Simulation Code/Methods/Tlasso.R")
                                    source("~/Simulation/Data Generation.R")
                                    Data <- Generate.Chi(n = n, dimen = dimen, Omega = Omega,
                                                         nu = nu, seed = seed*j, i.i.d = i.i.d,
                                                         variation = variation, nu.min = nu.min,
                                                         nu.max = nu.max)
                                    X <- Data$Cyclic.X;Vax <- Data$Cyclic.Vax
                                    K <- length(dimen);n.train <- dim(X)[K+1];n.val <- dim(Vax)[K+1]
                                    if(Heavy == "Matrix"){
                                      gamma <- array(0,c(dimen[1],n.train));gamma.vax <- array(0,c(dimen[1],n.val))
                                      for(j in 1:dimen[1]){
                                        for(i in 1:n.train){
                                          gamma[j,i] <- norm(X[j,,i],"2")/sqrt(dimen[2])
                                        }
                                        for(i in 1:n.val){
                                          gamma.vax[j,i] <- norm(Vax[j,,i],"2")/sqrt(dimen[2])
                                        }
                                      }
                                      x <- array(0,c(dimen,n.train));val <- array(0,c(dimen,n.val))
                                      for(i in 1:n.train){
                                        x[,,i] <- solve(diag(gamma[,i])) %*% X[,,i]
                                      }
                                      for(i in 1:n.val){
                                        val[,,i] <- solve(diag(gamma.vax[,i])) %*% Vax[,,i]
                                      }
                                    }else if(Heavy == "Constant"){
                                      gamma <- array(0,n.train);gamma.vax <- array(0,n.val)
                                      x <- array(0,c(dimen,n.train));val <- array(0,c(dimen,n.val))
                                      for(i in 1:n.train){
                                        gamma[i] <- norm(X[,,i],"F")/sqrt(prod(dimen))
                                        x[,,i] <- X[,,i] / gamma[i]
                                      }
                                      for(i in 1:n.val){
                                        gamma.vax[i] <- norm(Vax[,,i],"F")/sqrt(prod(dimen))
                                        val[,,i] <- Vax[,,i] / gamma.vax[i]
                                      }
                                    }else{
                                      x <- X;val <- Vax
                                    }
                                    if(penalty == 'adapt'){
                                      fit <- Tlasso.edit(x = x,val = val,T = t,lambda.min.ratio = rate[1],thres = thres)
                                      fit.adapt <- Tlasso.adapt(x = x,val = val, thres = thres,T = t,
                                                                Omega = fit$Omegahat,rate = rate[2])
                                      out <- Matrix.simulation.summary(Omega.hat.list = fit.adapt$Omegahat,
                                                                       Omega.true.list = Omega,
                                                                       norm.off.diag = norm.off.diag)
                                    }else{
                                      fit <- Tlasso.edit(x = x,val = val,T = t,lambda.min.ratio = rate[1])
                                      out <- Matrix.simulation.summary(Omega.hat.list = fit$Omegahat,
                                                                       Omega.true.list = Omega,
                                                                       norm.off.diag = norm.off.diag)
                                    }
                                    result <- list(Avg.F = out$av.error.f*100, Avg.Max = out$av.error.max*100,
                                                   Error.F = out$error.f*100, Error.Max = out$error.max*100,
                                                   Avg.Tpr = out$av.tpr*100, Avg.Fpr = out$av.fpr*100,
                                                   Avg.Tnr = out$av.tnr*100,
                                                   Tpr = out$tpr*100, Fpr = out$fpr*100,
                                                   Tnr = out$tnr * 100,
                                                   Optimality = fit$Optimality)
                                    result
                                  }
  Result <- as.data.frame(Cyclic)
  Error.F <- matrix(unlist(Result$Error.F),Run,byrow = TRUE)
  Error.Max <- matrix(unlist(Result$Error.Max),Run,byrow = TRUE)
  TPR <- matrix(unlist(Result$Tpr),Run, byrow = TRUE)
  FPR <- matrix(unlist(Result$Fpr),Run, byrow = TRUE)
  K = length(dimen)
  if(sum(sparse == seq(1:K))==K){
    Outcome <- list(Avg.Error.F.mean = round(mean(unlist(Result$Avg.F)),2), Avg.Error.F.sd =   round(sd(unlist(Result$Avg.F))/sqrt(Run),2),
                    Avg.Error.Max.mean = round(mean(unlist(Result$Avg.Max)),2), Avg.Error.Max.sd =   round(sd(unlist(Result$Avg.Max))/sqrt(Run),2),
                    Avg.TPR.mean = round(mean(unlist(Result$Avg.Tpr)),2), Avg.TPR.sd =   round(sd(unlist(Result$Avg.Tpr))/sqrt(Run),2),
                    Avg.FPR.mean = round(mean(unlist(Result$Avg.Fpr)),2), Avg.FPR.sd =   round(sd(unlist(Result$Avg.Fpr))/sqrt(Run),2),
                    Error.F.mean = round(colMeans(Error.F),2), Error.F.sd = round(apply(Error.F,2,sd)/sqrt(Run),2),
                    Error.Max.mean = round(colMeans(Error.Max),2), Error.Max.sd = round(apply(Error.Max,2,sd)/sqrt(Run),2),
                    TPR.mean = round(colMeans(TPR),2), TPR.sd = round(apply(TPR,2,sd)/sqrt(Run),2),
                    FPR.mean = round(colMeans(FPR),2), FPR.sd = round(apply(FPR,2,sd)/sqrt(Run),2))
  }
  if(sum(sparse == seq(1:K))!=K){
    Outcome <- list(Error.F.mean = round(colMeans(Error.F),2)[sparse], Error.F.sd = round(apply(Error.F,2,sd)/sqrt(Run),2)[sparse],
                    Error.Max.mean = round(colMeans(Error.Max),2)[sparse], Error.Max.sd = round(apply(Error.Max,2,sd)/sqrt(Run),2)[sparse],
                    TPR.mean = round(colMeans(TPR),2)[sparse], TPR.sd = round(apply(TPR,2,sd)/sqrt(Run),2)[sparse],
                    FPR.mean = round(colMeans(FPR),2)[sparse], FPR.sd = round(apply(FPR,2,sd)/sqrt(Run),2)[sparse])
  }
  return(Outcome)
}

Gemini.simul <- function(n, dimen, Omega,scale = NULL,
                         gam = 1,penalty = 'adapt', rate = c(1e-5,1e-6),
                         seed = 14, norm.off.diag = FALSE,
                         Heavy = "Matrix",nu = NULL,i.i.d = TRUE,
                         variation = TRUE, nu.min = NULL, nu.max = NULL,
                         Run = 100, sparse = c(1,2)){
  Gemini <- foreach(j = 1:Run, .combine =rbind,
                    .packages = c("tensr","ramify","stats",'rTensor','expm',
                                  "Tlasso","glasso","MASS"))%dopar%{
                                    source("~/Basic Setting.R")
                                    source("~/Simulation Evaluation.R")
                                    source("~/Simulation/Evaluation/Simulation Code/Methods/Gemini.R")
                                    source("~/Simulation/Data Generation.R")
                                    Data <- Generate.Chi(n = n, dimen = dimen, Omega = Omega,
                                                         nu = nu, seed = seed*j, i.i.d = i.i.d,
                                                         variation = variation, nu.min = nu.min,
                                                         nu.max = nu.max)
                                    X <- Data$Cyclic.X;Vax <- Data$Cyclic.Vax
                                    K <- length(dimen);n.train <- dim(X)[K+1];n.val <- dim(Vax)[K+1]
                                    if(Heavy == "Matrix"){
                                      gamma <- array(0,c(dimen[1],n.train));gamma.vax <- array(0,c(dimen[1],n.val))
                                      for(j in 1:dimen[1]){
                                        for(i in 1:n.train){
                                          gamma[j,i] <- norm(X[j,,i],"2")/sqrt(dimen[2])
                                        }
                                        for(i in 1:n.val){
                                          gamma.vax[j,i] <- norm(Vax[j,,i],"2")/sqrt(dimen[2])
                                        }
                                      }
                                      x <- array(0,c(dimen,n.train));val <- array(0,c(dimen,n.val))
                                      for(i in 1:n.train){
                                        x[,,i] <- solve(diag(gamma[,i])) %*% X[,,i]
                                      }
                                      for(i in 1:n.val){
                                        val[,,i] <- solve(diag(gamma.vax[,i])) %*% Vax[,,i]
                                      }
                                    }else if(Heavy == "Constant"){
                                      gamma <- array(0,n.train);gamma.vax <- array(0,n.val)
                                      x <- array(0,c(dimen,n.train));val <- array(0,c(dimen,n.val))
                                      for(i in 1:n.train){
                                        gamma[i] <- norm(X[,,i],"F")/sqrt(prod(dimen))
                                        x[,,i] <- X[,,i] / gamma[i]
                                      }
                                      for(i in 1:n.val){
                                        gamma.vax[i] <- norm(Vax[,,i],"F")/sqrt(prod(dimen))
                                        val[,,i] <- Vax[,,i] / gamma.vax[i]
                                      }
                                    }else{
                                      x <- X;val <- Vax
                                    }
                                    fit <- Gemini.partial(x = x,val = val,rate = rate[1])
                                    if(penalty == 'adapt'){
                                      fit.adapt <- Gemini.adapt(x = x,val = val, input = fit$corr,
                                                                inv.sd.row = fit$inv.sd.row,rate = rate[2])
                                      out <- Matrix.simulation.summary(Omega.hat.list = fit.adapt$Omegahat,
                                                                       Omega.true.list = Omega,
                                                                       norm.off.diag = norm.off.diag)
                                    }else{
                                      out <- Matrix.simulation.summary(Omega.hat.list = fit$Omegahat,
                                                                       Omega.true.list = Omega,
                                                                       norm.off.diag = norm.off.diag)
                                    }
                                    result <- list(Avg.F = out$av.error.f*100, Avg.Max = out$av.error.max*100,
                                                   Error.F = out$error.f*100, Error.Max = out$error.max*100,
                                                   Avg.Tpr = out$av.tpr*100, Avg.Fpr = out$av.fpr*100,
                                                   Avg.Tnr = out$av.tnr*100,
                                                   Tpr = out$tpr*100, Fpr = out$fpr*100,
                                                   Tnr = out$tnr * 100,
                                                   Optimality = fit$optimal)
                                    result
                                  }
  Result <- as.data.frame(Gemini)
  Error.F <- matrix(unlist(Result$Error.F),Run,byrow = TRUE)
  Error.Max <- matrix(unlist(Result$Error.Max),Run,byrow = TRUE)
  TPR <- matrix(unlist(Result$Tpr),Run, byrow = TRUE)
  FPR <- matrix(unlist(Result$Fpr),Run, byrow = TRUE)
  K = length(dimen)
  if(sum(sparse == seq(1:K))==K){
    Outcome <- list(Avg.Error.F.mean = round(mean(unlist(Result$Avg.F)),2), Avg.Error.F.sd =   round(sd(unlist(Result$Avg.F))/sqrt(Run),2),
                    Avg.Error.Max.mean = round(mean(unlist(Result$Avg.Max)),2), Avg.Error.Max.sd =   round(sd(unlist(Result$Avg.Max))/sqrt(Run),2),
                    Avg.TPR.mean = round(mean(unlist(Result$Avg.Tpr)),2), Avg.TPR.sd =   round(sd(unlist(Result$Avg.Tpr))/sqrt(Run),2),
                    Avg.FPR.mean = round(mean(unlist(Result$Avg.Fpr)),2), Avg.FPR.sd =   round(sd(unlist(Result$Avg.Fpr))/sqrt(Run),2),
                    Error.F.mean = round(colMeans(Error.F),2), Error.F.sd = round(apply(Error.F,2,sd)/sqrt(Run),2),
                    Error.Max.mean = round(colMeans(Error.Max),2), Error.Max.sd = round(apply(Error.Max,2,sd)/sqrt(Run),2),
                    TPR.mean = round(colMeans(TPR),2), TPR.sd = round(apply(TPR,2,sd)/sqrt(Run),2),
                    FPR.mean = round(colMeans(FPR),2), FPR.sd = round(apply(FPR,2,sd)/sqrt(Run),2))
  }
  if(sum(sparse == seq(1:K))!=K){
    Outcome <- list(Error.F.mean = round(colMeans(Error.F),2)[sparse], Error.F.sd = round(apply(Error.F,2,sd)/sqrt(Run),2)[sparse],
                    Error.Max.mean = round(colMeans(Error.Max),2)[sparse], Error.Max.sd = round(apply(Error.Max,2,sd)/sqrt(Run),2)[sparse],
                    TPR.mean = round(colMeans(TPR),2)[sparse], TPR.sd = round(apply(TPR,2,sd)/sqrt(Run),2)[sparse],
                    FPR.mean = round(colMeans(FPR),2)[sparse], FPR.sd = round(apply(FPR,2,sd)/sqrt(Run),2)[sparse])
  }
  return(Outcome)
}