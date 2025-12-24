library(MASS)
library(rTensor)
library(tensr)
library(expm)
library(stats)
source("~/Basic Setting.R")
source("~/Simulation/Hypothesis test/Distribution type/Normal.R")
source("~/Simulation/Hypothesis test/Distribution type/Contaminated.R")
source("~/Simulation/Hypothesis test/Distribution type/Chi.R")
source("~/Separate.R")

Test.Estimation <- function(n,dimen,alpha = 0.05,replicate = 500,
                            ncores = 96, Omega, 
                            seed = 14,nlambda  = c(20,10),
                            rate = c(1e-5,1e-6), norm.off.diag = F,
                            Simulation.type = c("Size","Power"),Test.type = c("Global","Local"),
                            alternative = c("two.sided","greater","less"),
                            dist.type = c("Normal","Chi","Contaminated"),
                            nu = NULL, i.i.d = NULL,variation = NULL, nu.min = NULL, nu.max = NULL,
                            value = NULL,prob = NULL,different = NULL){
  source("~/Basic Setting.R")
  source("~/Simulation/Hypothesis test/Distribution type/Normal.R")
  source("~/Simulation/Hypothesis test/Distribution type/Contaminated.R")
  source("~/Simulation/Hypothesis test/Distribution type/Chi.R")
  source("~/Separate.R")
  if(dist.type == "Normal"){
    hypothesis.result <- eval(parse(text = paste0(dist.type
                                                  ,'.Test(n = n,dimen = dimen,alpha = alpha,Omega = Omega,',
                                                  paste0('ncores = ncores,replicate = replicate'),
                                                  paste0(",seed = seed,nlambda = nlambda"),
                                                  paste0(",rate = rate, norm.off.diag = norm.off.diag"),
                                                  paste0(',Simulation.type = Simulation.type,Test.type = Test.type'),
                                                  paste0(',alternative = alternative)'))))
  }else if(dist.type == "Chi"){
    hypothesis.result <- eval(parse(text = paste0(dist.type
                                                  ,'.Test(n = n,dimen = dimen,alpha = alpha,Omega = Omega,',
                                                  paste0('ncores = ncores,replicate = replicate'),
                                                  paste0(",seed = seed,nlambda = nlambda"),
                                                  paste0(",rate = rate, norm.off.diag = norm.off.diag"),
                                                  paste0(',Simulation.type = Simulation.type,Test.type = Test.type'),
                                                  paste0(',nu = nu, i.i.d = i.i.d, variation = variation'),
                                                  paste0(',alternative = alternative'),
                                                  paste0(",nu.min = nu.min,nu.max = nu.max)"))))
  }else{
    hypothesis.result <- eval(parse(text = paste0(dist.type
                                                  ,'.Test(n = n,dimen = dimen,alpha = alpha,Omega = Omega,',
                                                  paste0('ncores = ncores,replicate = replicate'),
                                                  paste0(",seed = seed,nlambda = nlambda"),
                                                  paste0(",rate = rate, norm.off.diag = norm.off.diag"),
                                                  paste0(',Simulation.type = Simulation.type,Test.type = Test.type'),
                                                  paste0(',alternative= alternative'),
                                                  paste0(",value = value, prob = prob, different = different)"))))
  }
  return(hypothesis.result)
}

Test.Estimation.simul <- function(n,dimen,alpha = 0.05,replicate = 500,
                                  ncores = 96, Omega.list, Simul.index=seq(1:4),
                                  seed = 14,nlambda = c(20,10),
                                  rate = c(1e-5,1e-6), norm.off.diag = F,
                                  Simulation.type = c("Size","Power"),Test.type = c("Global","Local"),
                                  alternative = c("two.sided","greater","less"),
                                  dist.type = c("Normal","Chi","Contaminated"),
                                  nu = NULL, i.i.d = NULL,variation = NULL, nu.min = NULL, nu.max = NULL,
                                  value = NULL,prob = NULL,different = NULL){
  q <- nrow(dimen)
  if(dist.type == "Normal" | Simulation.type == "Size"){
    Size <- array(0,q)
    for(i in 1:q){
      dim1 <- dimen[i,]
      if(dist.type == "Normal"){
        Size[i] <- eval(parse(text = paste0('Test.Estimation(n = n,dimen = dim1,alpha = alpha,Omega = Omega.list[[i]],',
                                            paste0('ncores = ncores,replicate = replicate'),
                                            paste0(",seed = seed,nlambda = nlambda"),
                                            paste0(",rate = rate, norm.off.diag = norm.off.diag"),
                                            paste0(',Simulation.type = Simulation.type,Test.type = Test.type'),
                                            paste0(',alternative = alternative, dist.type = dist.type,'),
                                            paste0("nu = nu, i.i.d = i.i.d, variation = variation,"),
                                            paste0("nu.min = nu.min, nu.max = nu.max,"),
                                            paste0("value = value, prob = prob, different = different)"))))
      }else{
        A <- eval(parse(text = paste0('Test.Estimation(n = n,dimen = dim1,alpha = alpha,Omega = Omega.list[[i]],',
                                      paste0('ncores = ncores,replicate = replicate'),
                                      paste0(",seed = seed,nlambda = nlambda"),
                                      paste0(",rate = rate, norm.off.diag = norm.off.diag"),
                                      paste0(',Simulation.type = Simulation.type,Test.type = Test.type'),
                                      paste0(',alternative = alternative, dist.type = dist.type,'),
                                      paste0("nu = nu, i.i.d = i.i.d, variation = variation,"),
                                      paste0("nu.min = nu.min, nu.max = nu.max,"),
                                      paste0("value = value, prob = prob, different = different)"))))
        Size[i] <- A$Size
      } 
    }
    Power <- NULL
  }else{
    Power <- array(0,q)
    for(i in 1:q){
      dim1 <- dimen[i,]
      A <- eval(parse(text = paste0('Test.Estimation(n = n,dimen = dim1,alpha = alpha,Omega = Omega.list[[i]],',
                                    paste0('ncores = ncores,replicate = replicate'),
                                    paste0(",seed = seed,nlambda = nlambda"),
                                    paste0(",rate = rate, norm.off.diag = norm.off.diag"),
                                    paste0(',Simulation.type = Simulation.type,Test.type = Test.type'),
                                    paste0(',alternative = alternative, dist.type = dist.type,'),
                                    paste0("nu = nu, i.i.d = i.i.d, variation = variation,"),
                                    paste0("nu.min = nu.min, nu.max = nu.max,"),
                                    paste0("value = value, prob = prob, different = different)"))))
      Power[i] <- A$Power
    }
    Size <- NULL
  }
  if(q == 8){
    Simul.index <- seq(1:q)
  }
  if(Simulation.type == "Size"){
    names(Size) = c("C1","C2","R1","R2")[Simul.index]
  }else{
    names(Power) = c("C1","C2","R1","R2")[Simul.index]
  }
  if(Test.type == "Global"){
    Test <- paste0("H0: The Data matrix follows MGGM. vs ",
                   paste0("H1: Not all the heavy-tails are constant."))
  }else{
    Test <- paste0("H0: The Data matrix follows U-MANGO. vs ",
                   paste0("H1: The Data matrix follows H-MANGO."))
  }
  if(dist.type == "Normal"){
    Distribution <- "Matrix Normal Distribution"
  }else if(dist.type == "Chi"){
    Distribution <- "Chi-Squared distribution."
    if(isFALSE(variation)){
      Distribution <- paste0(paste0(Distribution, " Constant heavy-tail. "),
                             paste0("Degree of freedom is ",paste0(nu,".")))
    }else{
      if(i.i.d){
        Distribution <- paste0(paste0(Distribution, " Row heavy-tail. "),
                               paste0("Degree of freedom is ",paste0(nu,".")))
      }else{
        Distribution <- paste0(paste0(Distribution, " Row heavy-tail. "),
                               paste0("The heavy-tail doesn't follow i.i.d. assumption. "),
                               paste0("The degree of freedom follows Unif[",nu.min,",",nu.max,"]."))
      }
    }
  }else{
    Distribution <- "Contaminated Random Variable,"
    if(different){
      Distribution <- paste0(paste0(Distribution," Row heavy-tail."),
                             paste0(" Heavy tail is generated by discrete random variable,"),
                             paste0(" where ",paste0("P(gamma=",value, ")=",prob,collapse = ", ")),'.')
    }else{
      Distribution <- paste0(paste0(Distribution," Constant heavy-tail."),
                             paste0(" Heavy tail is generated by discrete random variable,"),
                             paste0(" where ",paste0("P(gamma=",value, ")=",prob,collapse = ", ")),'.')
    }
  }
  result <- list(Size = Size, Power = Power, 
                 Distribution = Distribution,Test = Test,
                 alternative = alternative,
                 Simulation.type = Simulation.type)
  return(result)
}