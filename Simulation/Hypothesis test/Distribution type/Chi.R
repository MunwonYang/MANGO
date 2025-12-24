source("~/Basic Setting.R")
source("~/Separate.R")
source("~/helper.R")
library(doParallel)
library(foreach)
library(MASS)
library(rTensor)
library(tensr)
library(expm)
library(stats)

Chi.Test <- function(n,dimen,alpha = 0.05, replicate = 1000,nlambda = c(50,20),
                     ncores = 96, seed = 14, Omega,
                     norm.off.diag = FALSE, rate = c(1e-5,1e-6),
                     Simulation.type = c("Size","Power"),
                     alternative = c("two.sided","greater","less"),
                     Test.type = c("Global",'Local'),
                     nu = NULL, i.i.d = NULL, variation = NULL,
                     nu.min = NULL,nu.max = NULL){
  dist.type <- "Chi"
  if(Test.type == "Global"){
    Simulation.type <- "Power"
    if(length(alternative)!=1){
      alternative <- "greater"
    }
    Test.stat <- foreach(iter = 1:replicate,.combine = rbind,
                         .packages = c("tensr","ramify","stats",'rTensor','expm',
                                       "Tlasso","glasso","MASS"))%dopar%{
                                         source("~/Basic Setting.R")
                                         source("~/Separate.R")
                                         source("~/helper.R")
                                         seed.value<- seed* iter
                                         Data<- Generate.Chi(n = n,dimen = dimen,Omega = Omega,seed = seed.value,
                                                              nu = nu, variation = variation, i.i.d = i.i.d,
                                                              nu.min = nu.min, nu.max = nu.max)
                                         X <- Data$x;K <- length(dimen)
                                         fit <- cv.Separate(x = X, seed = seed.value,
                                                            nlambda = nlambda,rate = rate[1])
                                         fit.adapt <- cv.Separate.penalty(x = X,seed = seed.value,
                                                                          Omega = fit$Omegahat,
                                                                          nlambda = nlambda,penalty = "adapt",
                                                                          rate = rate[2])
                                         w <- atrans(X,lapply(fit.adapt$Omegahat, function(x) sqrtm(x)))
                                         pval <- Matrix.power(x = w, type = Test.type, alpha = alpha)
                                         result <- list(statistic = pval$statistic,
                                                        power = pval$power)
                                         return(result)
                                       }
    power <- as.data.frame(Test.stat);Power = mean(unlist(power$power));Size <- NULL
  }else{
    eval(parse(text = paste0("T.star=Hypothesis.test.parameter(n = n,dimen = dimen,",
                             "ncores = ncores)")))
    H0.mu <- T.star$mu;H0.sigma <- T.star$sigma
    if(length(alternative)!=1){
      alternative <- "two.sided"
    }
    if(Simulation.type == "Size"){
      Test.stat <- foreach(iter = 1:replicate,.combine = rbind,
                           .packages = c("tensr","ramify","stats",'rTensor','expm',
                                         "Tlasso","glasso","MASS"))%dopar%{
                                           source("~/Basic Setting.R")
                                           source("~/Separate.R")
                                           source("~/helper.R")
                                           seed.value<- seed* iter
                                           Data<- Generate.Chi(n = n,dimen = dimen,Omega = Omega,seed = seed.value,
                                                               nu = nu, variation = variation, i.i.d = i.i.d,
                                                               nu.min = nu.min, nu.max = nu.max)
                                           X <- Data$x;K <- length(dimen);con.test<- Matrix.Surrogate(X,type = "Constant");W <- sweep(X,3,con.test,"/")
                                           fit <- cv.Separate(x = W, seed = seed.value,
                                                              nlambda = nlambda,rate = rate[1])
                                           fit.adapt <- cv.Separate.penalty(x = W,seed = seed.value,
                                                                            Omega = fit$Omegahat,
                                                                            nlambda = nlambda,penalty = "adapt",
                                                                            rate = rate[2])
                                           w <- atrans(W,lapply(fit.adapt$Omegahat, function(x) sqrtm(x)))
                                           mat <- Matrix.Surrogate(w,type = "Matrix");con <- Matrix.Surrogate(w,type = "Constant")
                                           est <- mean((sweep(mat^2,2,con^2,"/")-1)^2);test.stat = (est - H0.mu) / H0.sigma
                                           if(alternative == "two.sided"){
                                             if(test.stat <0){
                                               p.value = pnorm(test.stat) * 2
                                             }else{
                                               p.value = pnorm(test.stat,lower.tail = F) * 2
                                             }
                                           }else{
                                             if(test.stat <0){
                                               p.value = pnorm(test.stat)
                                             }else{
                                               p.value = pnorm(test.stat,lower.tail = F)
                                             }
                                           }
                                           return(p.value)
                                         }
      Size <- sum(Test.stat <= alpha) / replicate;Power <- NULL
    }else{
      Test.stat <- foreach(iter = 1:replicate,.combine = rbind,
                           .packages = c("tensr","ramify","stats",'rTensor','expm',
                                         "Tlasso","glasso","MASS"))%dopar%{
                                           source("~/Basic Setting.R")
                                           source("~/Separate.R")
                                           source("~/helper.R")
                                           seed.value<- seed* iter
                                           Data<- Generate.Chi(n = n,dimen = dimen,Omega = Omega,seed = seed.value,
                                                               nu = nu, variation = variation, i.i.d = i.i.d,
                                                               nu.min = nu.min, nu.max = nu.max)
                                           X <- Data$x;K <- length(dimen);con.test<- Matrix.Surrogate(X,type = "Constant");W <- sweep(X,3,con.test,"/")
                                           fit <- cv.Separate(x = W, seed = seed.value,
                                                              nlambda = nlambda,rate = rate[1])
                                           fit.adapt <- cv.Separate.penalty(x = W,seed = seed.value,
                                                                            Omega = fit$Omegahat,
                                                                            nlambda = nlambda,penalty = "adapt",
                                                                            rate = rate[2])
                                           w <- atrans(W,lapply(fit.adapt$Omegahat, function(x) sqrtm(x)))
                                           mat <- Matrix.Surrogate(w,type = "Matrix");con <- Matrix.Surrogate(w,type = "Constant")
                                           est <- mean((sweep(mat^2,2,con^2,"/")-1)^2);test.stat = (est - H0.mu) / H0.sigma
                                           if(alternative == "two.sided"){
                                             power <- pnorm(qnorm(alpha/2) + abs(test.stat))
                                           }else{
                                             if(alternative == "greater"){
                                               power <- pnorm(qnorm(alpha) + test.stat)
                                             }else{
                                               power <- 1 - pnorm(qnorm(1-alpha) + test.stat)
                                             }
                                           }
                                           return(power)
                                         }
      Power <- mean(Test.stat);Size <- NULL
    }
  }
  result <- list(Size = Size, Power = Power)
  return(result)
}