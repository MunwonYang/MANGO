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

Normal.Test <- function(n,dimen,alpha = 0.05, replicate = 1000,nlambda = c(50,20),
                        ncores = 96, seed = 14, Omega,
                        norm.off.diag = FALSE, rate = c(1e-5,1e-6),
                        Simulation.type = "Size",
                        Test.type = c("Global","Local"),
                        alternative = c("two.sided","greater","less")){
  dist.type <- "Normal"
  Simulation.type <- "Size"
  if(Test.type == "Global"){
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
                                         Data <- Generate.Normal(n = n,dimen = dimen, Omega = Omega, seed = seed.value)
                                         X <- Data$x;K <- length(dimen)
                                         fit <- cv.Separate(x = X, seed = seed.value,
                                                            nlambda = nlambda,rate = rate[1])
                                         fit.adapt <- cv.Separate.penalty(x = X,seed = seed.value,
                                                                          Omega = fit$Omegahat,
                                                                          nlambda = nlambda,penalty = "adapt",
                                                                          rate = rate[2])
                                         w <- atrans(X,lapply(fit.adapt$Omegahat, function(x) sqrtm(x)))
                                         pval <- Matrix.tail.test(w,type = Test.type, alternative = alternative)
                                         result <- list(statistic = pval$statistic,
                                                        pvalue = pval$p.value)
                                         return(result)
                                       }
    size <- as.data.frame(Test.stat)
    Size <- sum(size$pvalue<=alpha) / replicate
  }else{
    eval(parse(text = paste0("T.star=Hypothesis.test.parameter(n = n,dimen = dimen,",
                             "ncores = ncores)")))
    if(length(alternative)!=1){
      alternative <- "two.sided"
    }
    H0.mu <- T.star$mu;H0.sigma <- T.star$sigma
    Test.stat <- foreach(iter = 1:replicate,.combine = rbind,
                         .packages = c("tensr","ramify","stats",'rTensor','expm',
                                       "Tlasso","glasso","MASS"))%dopar%{
                                         source("~/Basic Setting.R")
                                         source("~/Separate.R")
                                         source("~/helper.R")
                                         seed.value<- seed* iter
                                         Data <- Generate.Normal(n = n,dimen = dimen, Omega = Omega, seed = seed.value)
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
    Size <- sum(Test.stat <= alpha) / replicate
  }
  return(Size)
}


